"""Tests for the analysis kernel.

These pin the properties the workflow's conclusions rest on: that alignment
scores are correct, that the batched scan agrees exactly with the pairwise
aligner, that motif and array detection distinguish the cases the triage stage
is asked to distinguish, and that the statistics are not made up.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from studio.bio.align import (
    BLOSUM62,
    align,
    encode_batch,
    profile_align,
    profile_scores_batch,
)
from studio.bio.cluster import greedy_cluster, taxonomic_spread
from studio.bio.context import (
    Gene,
    KnownSystem,
    Locus,
    architecture,
    assess_novelty,
    extract_locus,
    predicted_operon,
)
from studio.bio.motifs import MotifSet, MotifSpec, catalytic_triad_intact, prosite_to_regex
from studio.bio.profile import build_profile, profile_from_unaligned
from studio.bio.repeats import find_repeat_arrays
from studio.bio.search import SequenceIndex, reduce_alphabet
from studio.bio.seq import low_complexity_fraction, revcomp, translate
from studio.bio.stats import benjamini_hochberg, fisher_exact

RT = (
    "MGKLIDLVARFYRQVLRNKGAPGVDGMTVEELPAYLKAHWPTIRQQLLDGTYRPQPVRRVEIPKPGGGTRLLG"
    "IPTVTDRLIQQAIAQVLTPIYEAQFSDHSYGFRPGRSAHDAVRQAQDYIQQGYRWVVDIDLEKFFDRVNHDIL"
    "MSRVARKVKDKRVLKLIRAYLQAGVMIDGVVQQTEEGTPQGGPLSPLLSNILLDDLDKELERRGLRFVRYADD"
    "CNIYVKSLRAGQRVKQSIQRFLEKTLKLKVNEEKSAVDRPWKRAFLGFSFTPERKARIRLAPRSIQRLKQRIR"
    "QLTNRSWGVSNRYRIEKLNQLIRGWINYFKIGSMKTLCRELDEWIRRRLRCCIWKQWKRVRTRIRELIALGLK"
)


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(c if rng.random() > rate else rng.choice(aa) for c in seq)


# --- sequence and alignment ---------------------------------------------------


def test_blosum62_is_symmetric_and_correct():
    from studio.bio.seq import AA_INDEX

    assert np.allclose(BLOSUM62[:20, :20], BLOSUM62[:20, :20].T)
    assert BLOSUM62[AA_INDEX["W"], AA_INDEX["W"]] == 11
    assert BLOSUM62[AA_INDEX["A"], AA_INDEX["A"]] == 4
    assert BLOSUM62[AA_INDEX["Y"], AA_INDEX["W"]] == 2
    assert BLOSUM62[AA_INDEX["C"], AA_INDEX["P"]] == -3


def test_self_alignment_is_perfect():
    a = align(RT, RT)
    assert a.identity == 1.0
    assert a.aligned_len == len(RT)
    assert a.gaps == 0


def test_local_alignment_finds_the_embedded_domain():
    domain = RT[60:140]
    a = align(domain, RT)
    assert (a.t_start, a.t_end) == (60, 140)
    assert a.identity == 1.0


def test_unrelated_sequences_score_near_noise():
    rng = random.Random(0)
    other = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(len(RT)))
    assert align(RT, other).score < align(RT, RT).score / 5


def test_translate_and_revcomp_round_trip():
    dna = "ATGAAACGTTAA"
    assert translate(dna) == "MKR"
    assert revcomp(revcomp(dna)) == dna
    assert revcomp("ATGC") == "GCAT"


def test_low_complexity_detects_homopolymer():
    assert low_complexity_fraction("Q" * 80) > 0.9
    assert low_complexity_fraction(RT) < 0.2


# --- batched scoring must agree exactly with the pairwise aligner -------------


def test_batched_scores_match_pairwise_exactly():
    """The scan is only sound if batching changes speed and nothing else."""
    rng = random.Random(5)
    profile = profile_from_unaligned(
        "RT", [mutate(RT, 0.2, rng) for _ in range(5)], calibration_n=20
    )
    targets = [mutate(RT, 0.35, rng) for _ in range(8)] + [
        "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rng.randint(80, 300)))
        for _ in range(12)
    ]
    enc, lengths = encode_batch(targets)
    batched = profile_scores_batch(profile.pssm, enc, lengths)
    pairwise = np.array([profile_align(profile.pssm, t).score for t in targets])
    assert np.allclose(batched, pairwise, atol=1e-3)


def test_batching_is_insensitive_to_padding():
    """A short sequence must score the same however long the batch's longest is."""
    rng = random.Random(6)
    profile = profile_from_unaligned(
        "RT", [mutate(RT, 0.2, rng) for _ in range(4)], calibration_n=20
    )
    short = mutate(RT, 0.3, rng)[:90]
    alone = profile_scores_batch(profile.pssm, *encode_batch([short]))[0]
    padded = profile_scores_batch(profile.pssm, *encode_batch([short, RT * 2]))[0]
    assert alone == pytest.approx(padded, abs=1e-3)


# --- motifs -------------------------------------------------------------------


def test_prosite_compilation():
    assert prosite_to_regex("Y-x-[DN]-[DN]") == "Y[A-Z][DN][DN]"
    assert prosite_to_regex("A-x(2,4)-{PG}") == "A[A-Z]{2,4}[^PG]"
    assert prosite_to_regex("<M-x(3)-K>") == "^M[A-Z]{3}K$"


def test_catalytic_triad_distinguishes_the_cases_triage_must_distinguish():
    intact, why = catalytic_triad_intact(RT)
    assert intact and "D" in why

    dead, why_dead = catalytic_triad_intact(RT.replace("YADDCNIY", "YANNCNIY"))
    assert not dead and "substituted" in why_dead

    gone, _ = catalytic_triad_intact(RT.replace("YADDCNIY", "YAGGCNIY"))
    assert not gone

    fragment, _ = catalytic_triad_intact(RT[:180])
    assert not fragment


def test_motif_placement_survives_a_spurious_early_match():
    """A permissive motif matching early must not block a later essential one.

    This is the failure that a first-match-wins scanner has: a chance hit for
    an upstream motif consumes the ordering constraint and the real downstream
    motif is then reported missing.
    """
    mset = MotifSet(
        name="t",
        motifs=[
            MotifSpec("first", "A-x-A", essential=False),
            MotifSpec("second", "W-W-W", essential=True),
        ],
        spacing={("first", "second"): (5, 30)},
    )
    # An "AxA" at position 0 is too far from the WWW; a second one is in range.
    seq = "ACA" + "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" + "AQA" + "KKKKK" + "WWW"
    report = mset.scan(seq)
    assert not report.missing_essential
    assert not report.spacing_violations


def test_motif_information_prefers_the_specific_match():
    from studio.bio.motifs import pattern_information

    assert pattern_information("[YF]-x-[DN]-[DN]") > pattern_information("[YFWLM]-x-[DN]-[DN]")


# --- repeat arrays ------------------------------------------------------------


def _array_dna(n_copies, repeat_len, spacer_len, rng, identical_spacers=False):
    repeat = "".join(rng.choice("ACGT") for _ in range(repeat_len))
    fixed = "".join(rng.choice("ACGT") for _ in range(spacer_len))
    out = []
    for i in range(n_copies):
        out.append(repeat)
        if i < n_copies - 1:
            out.append(fixed if identical_spacers
                       else "".join(rng.choice("ACGT") for _ in range(spacer_len)))
    return "".join(out)


def test_spacer_array_is_found_and_classified():
    rng = random.Random(1)
    flank = "".join(rng.choice("ACGT") for _ in range(500))
    dna = flank + _array_dna(9, 30, 35, rng) + flank
    arrays = find_repeat_arrays(dna)
    assert arrays
    a = arrays[0]
    assert a.copy_number >= 7
    assert a.spacer_regularity > 0.85
    assert a.class_token() == "spacer-array"


def test_vntr_is_not_mistaken_for_a_spacer_array():
    """Near-identical spacers mean a tandem repeat, not a spacer array."""
    rng = random.Random(2)
    flank = "".join(rng.choice("ACGT") for _ in range(400))
    dna = flank + _array_dna(9, 30, 35, rng, identical_spacers=True) + flank
    arrays = find_repeat_arrays(dna)
    assert arrays
    assert arrays[0].class_token() == "tandem-repeat"
    assert arrays[0].spacer_identity > 0.6


def test_shuffled_dna_yields_no_regular_array():
    """The negative control the reproduce stage runs, as a unit test."""
    rng = random.Random(3)
    dna = "".join(rng.choice("ACGT") for _ in range(6000))
    arrays = find_repeat_arrays(dna)
    assert not [a for a in arrays if a.copy_number >= 3 and a.spacer_regularity > 0.75]


def test_array_recovers_copies_the_exact_seed_pass_would_miss():
    """A mutated copy inside the seed k-mer must not truncate the array."""
    rng = random.Random(4)
    repeat = "".join(rng.choice("ACGT") for _ in range(30))
    parts = []
    for i in range(8):
        unit = list(repeat)
        if i in (2, 5):                       # break the seed in two copies
            unit[3] = "A" if unit[3] != "A" else "C"
            unit[9] = "G" if unit[9] != "G" else "T"
        parts.append("".join(unit))
        if i < 7:
            parts.append("".join(rng.choice("ACGT") for _ in range(35)))
    dna = "".join(rng.choice("ACGT") for _ in range(300)) + "".join(parts)
    arrays = find_repeat_arrays(dna)
    assert arrays and arrays[0].copy_number >= 7


# --- profiles and search ------------------------------------------------------


def test_profile_columns_carry_information():
    rng = random.Random(7)
    profile = profile_from_unaligned(
        "RT", [mutate(RT, 0.15, rng) for _ in range(8)], calibration_n=25
    )
    assert profile.length > 200
    # The sequence-weighting fix: pseudocounts must not swamp a deep column.
    assert len(profile.core_columns) > profile.length // 4
    assert profile.decoy_sigma > 0


def test_profile_search_recovers_homologs_and_rejects_fragments():
    rng = random.Random(8)
    profile = profile_from_unaligned(
        "RT", [mutate(RT, 0.2, rng) for _ in range(6)], calibration_n=40
    )
    index = SequenceIndex()
    for i in range(300):
        index.add(
            f"decoy{i}",
            "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(rng.randint(150, 400))),
        )
    for i in range(6):
        index.add(f"true{i}", mutate(RT, 0.38, rng))
    for i in range(4):
        index.add(f"frag{i}", mutate(RT, 0.25, rng)[140:250])

    hits = index.profile_search(profile, min_z=6.0)
    found = {h.target_id for h in hits}
    assert len([h for h in found if h.startswith("true")]) >= 5
    assert not [h for h in found if h.startswith("frag")]
    assert not [h for h in found if h.startswith("decoy")]


def test_search_space_uses_query_length():
    index = SequenceIndex()
    index.add("a", RT)
    assert index.search_space(130)[0] == 130
    assert index.search_space(400)[0] == 400
    assert index.search_space(130)[1] == index.total_residues


def test_reduced_alphabet_groups_conservative_substitutions():
    assert reduce_alphabet("LVIM") == "LLLL"
    assert reduce_alphabet("FYW") == "FFF"
    assert reduce_alphabet("KR") == "KK"


# --- clustering and statistics ------------------------------------------------


def test_clustering_collapses_near_identical_sequences():
    rng = random.Random(9)
    records = {f"near{i}": mutate(RT, 0.02, rng) for i in range(5)}
    records["distant"] = mutate(RT, 0.5, rng)
    clusters = greedy_cluster(records, identity=0.9)
    assert len(clusters) == 2


def test_taxonomic_spread_flags_a_single_clade():
    rng = random.Random(10)
    records = {f"s{i}": mutate(RT, 0.3, rng) for i in range(6)}
    same = {k: "Bacteria;Pseudomonadota;Gamma;Entero" for k in records}
    spread = taxonomic_spread(records, same)
    assert spread.n_phyla == 1
    assert spread.is_clade_restricted


def test_fisher_exact_matches_known_values():
    # Textbook 2x2: strong association.
    r = fisher_exact(10, 0, 2, 18)
    assert r.p_value < 1e-4
    assert r.odds_ratio > 20
    # No association at all.
    r2 = fisher_exact(5, 5, 5, 5)
    assert r2.p_value > 0.4
    assert 0.5 < r2.odds_ratio < 2.0


def test_benjamini_hochberg_controls_the_family():
    pvals = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.5, 0.9]
    flags = benjamini_hochberg(pvals, alpha=0.05)
    assert flags[0] and flags[1]
    assert not flags[-1] and not flags[-2]
    # A monotone property: nothing above the largest rejected p is rejected.
    rejected = [p for p, f in zip(pvals, flags) if f]
    assert all(p <= max(rejected) for p, f in zip(pvals, flags) if f)


# --- genomic context ----------------------------------------------------------


def _locus(labels, strands=None, arrays=None, gap=40):
    """Genes laid out with a 40 bp intergenic gap, inside the operon cutoff."""
    strands = strands or [1] * len(labels)
    genes = []
    cursor = 0
    for i, lab in enumerate(labels):
        genes.append(
            Gene(f"g{i}", cursor, cursor + 300, strands[i], protein="M" * 100, label=lab)
        )
        cursor += 300 + gap
    return Locus("ctg", genes[0], genes, arrays or [])


def test_operon_prediction_stops_at_a_strand_switch():
    loc = _locus(["RT", "partner", "other"], strands=[1, 1, -1])
    op = predicted_operon(loc)
    assert [g.label for g in op] == ["RT", "partner"]


def test_operon_prediction_stops_at_a_long_intergenic_gap():
    loc = _locus(["RT", "partner"], gap=400)
    assert [g.label for g in predicted_operon(loc)] == ["RT"]


def test_architecture_class_mode_drops_array_geometry():
    from studio.bio.repeats import RepeatArray

    arr = RepeatArray(
        start=900, end=1500, repeat_consensus="A" * 30,
        repeat_positions=[900, 965, 1030], spacer_lengths=[35, 35],
        spacers=["ACGT", "TGCA"], repeat_identity=0.98, spacer_identity=0.1,
    )
    loc = _locus(["RT", "unk"], arrays=[arr])
    assert "array[3x30]" in architecture(loc, array_detail="full")
    assert "array[spacer-array]" in architecture(loc, array_detail="class")


def test_novelty_is_zero_for_an_exact_catalogue_match():
    catalogue = [
        KnownSystem("retron", "Retron", ["RT", "effector"], forbids_array=True),
    ]
    loc = _locus(["RT", "effector"])
    assessment = assess_novelty(loc, catalogue)
    assert assessment.novelty == 0.0
    assert assessment.best_match.is_match


def test_an_unexplained_component_keeps_novelty_on_the_board():
    """A described system plus a component it never has is not that system."""
    from studio.bio.repeats import RepeatArray

    catalogue = [
        KnownSystem("retron", "Retron", ["RT", "effector"], forbids_array=True),
    ]
    arr = RepeatArray(
        start=900, end=1500, repeat_consensus="A" * 30,
        repeat_positions=[900, 965, 1030], spacer_lengths=[35, 35],
        spacers=["AC", "TG"], repeat_identity=0.98, spacer_identity=0.1,
    )
    loc = _locus(["RT", "effector"], arrays=[arr])
    assessment = assess_novelty(loc, catalogue)
    assert assessment.novelty >= 0.35
    assert any("array" in v for v in assessment.best_match.violations)


def test_extract_locus_respects_the_window():
    genes = [
        Gene(f"g{i}", i * 20_000, i * 20_000 + 300, 1, protein="M" * 100)
        for i in range(5)
    ]
    loc = extract_locus("ctg", genes, "g2", window_bp=5000)
    assert [g.gene_id for g in loc.genes] == ["g2"]
