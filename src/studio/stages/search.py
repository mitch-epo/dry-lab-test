"""Stage 3 -- search: find family members that fit no described system.

"It then searches for family members or genomic neighbors that fit no
described system..."

The stage does four things in order, and the order matters:

1. scan the database with the calibrated profile, keeping the funnel's
   rejection counts so a short list can be shown to be the product of
   filtering rather than of a broken index;
2. pull each hit's genomic neighbourhood and label the neighbours by homology
   to the reference set, leaving genuinely unmatched genes as ``unk``;
3. score each locus against every described system, so novelty is a measured
   quantity with a named runner-up rather than an impression; and
4. group loci by architecture and compute, per group, the statistics that
   decide whether an association is real -- partner conservation across
   independent loci, enrichment against a background drawn from the same
   genomes, sequence redundancy, and taxonomic spread.

Step 4 is what separates this from a list of interesting-looking loci. A single
locus can never support a claim about a system; a recurrent architecture can.
"""

from __future__ import annotations

import time
from collections import defaultdict
from dataclasses import asdict, dataclass, field

from ..bio.cluster import taxonomic_spread
from ..bio.context import (
    Locus,
    architecture,
    assess_novelty,
    conservation_across_loci,
    core_architecture,
    predicted_operon,
)
from ..bio.motifs import catalytic_triad_intact
from ..bio.profile import Profile
from ..bio.search import Hit, SequenceIndex
from ..bio.seq import low_complexity_fraction
from ..bio.stats import benjamini_hochberg, neighbour_enrichment
from ..campaign import CampaignConfig, Run
from ..evidence import strength_for
from ..schemas import Candidate, Evidence, StageRecord, Strength
from ..sources.base import DataSource
from .survey import FamilyLabeller, build_reference_index


@dataclass
class ArchitectureGroup:
    """All loci sharing one core architecture -- the unit a claim is made about."""

    signature: str
    members: list[str] = field(default_factory=list)
    partner_labels: list[str] = field(default_factory=list)
    conservation: dict[str, float] = field(default_factory=dict)
    n_clusters: int = 0
    n_phyla: int = 0
    clade_restricted: bool = True
    enrichment: dict[str, dict] = field(default_factory=dict)
    spread_summary: str = ""


@dataclass
class SearchResult:
    scanned: int = 0
    funnel: str = ""
    profile_hits: int = 0
    loci: int = 0
    novel_loci: int = 0
    groups: list[dict] = field(default_factory=list)
    candidates: int = 0
    rejected: dict[str, int] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


def _array_quality(locus: Locus) -> tuple[float, str, object]:
    """Score the best repeat array at a locus, and say why.

    A well-formed spacer array -- three or more copies, regular spacing,
    distinct spacers -- is a strong architectural signal. A tandem repeat with
    near-identical spacers is a VNTR and carries none of that meaning, so it
    scores near zero rather than being silently counted as an array.
    """
    best, best_q, why = None, 0.0, "no repeat array in the neighbourhood"
    for a in locus.arrays:
        if a.copy_number < 3:
            continue
        distinctness = max(0.0, 1.0 - a.spacer_identity / 0.6)
        q = (
            0.45 * min(1.0, a.copy_number / 8.0)
            + 0.30 * a.spacer_regularity
            + 0.25 * distinctness
        ) * (0.15 if a.spacer_identity > 0.6 else 1.0)
        if q > best_q:
            best, best_q = a, q
            why = a.summary() + f" -- {a.classify()}"
    return best_q, why, best


def _operon_distance(locus: Locus, array) -> int:
    """Base pairs between the predicted operon and an array."""
    if array is None:
        return 10 ** 9
    op = predicted_operon(locus)
    if not op:
        return 10 ** 9
    lo, hi = min(g.start for g in op), max(g.end for g in op)
    if array.start > hi:
        return array.start - hi
    if array.end < lo:
        return lo - array.end
    return 0


def run_search(
    run: Run,
    config: CampaignConfig,
    source: DataSource,
    profile: Profile,
    index: SequenceIndex,
    systems: list,
    labeller: FamilyLabeller,
) -> tuple[SearchResult, list[Candidate]]:
    """Execute the search. Returns the summary and the candidate list."""
    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    rejected: dict[str, int] = defaultdict(int)

    # -- 1: profile scan ----------------------------------------------------
    hits: list[Hit] = index.profile_search(
        profile,
        max_evalue=config.profile_max_evalue,
        min_z=config.profile_min_z,
        min_core_coverage=config.profile_min_core_coverage,
    )
    funnel = index.last_scan.summary()

    # -- 2: loci and labels -------------------------------------------------
    label_cache: dict = {}
    loci: dict[str, Locus] = {}
    hit_by_id: dict[str, Hit] = {h.target_id: h for h in hits}
    for h in hits:
        loc = source.locus(h.target_id, window_bp=config.neighbourhood_bp)
        if loc is None:
            rejected["no genomic context available"] += 1
            continue
        calls = labeller.label_genes(
            loc.genes,
            min_z=config.label_min_z,
            min_coverage=config.label_min_coverage,
            cache=label_cache,
        )
        for g in loc.genes:
            call = calls.get(g.gene_id)
            g.label = call.label if call else "unk"
            if call and call.label != "unk":
                g.annotation = g.annotation or f"{call.label} (z={call.z:.0f})"
        loc.anchor.label = config.family
        loci[h.target_id] = loc

    # -- 3: novelty against the catalogue -----------------------------------
    assessments = {pid: assess_novelty(loc, systems) for pid, loc in loci.items()}
    novel = {pid: a for pid, a in assessments.items() if a.novelty >= config.min_novelty}
    for pid, a in assessments.items():
        if pid not in novel:
            rejected[f"explained by {a.best_match.name if a.best_match else 'catalogue'}"] += 1

    # -- 4: group by architecture and compute group statistics --------------
    # Group on the class-level signature, not the geometry: see the note in
    # studio.bio.context.architecture. The full-detail string is kept on each
    # candidate for the report.
    groups: dict[str, ArchitectureGroup] = {}
    for pid in novel:
        sig = core_architecture(
            loci[pid], max_gap=config.operon_max_gap, array_detail="class"
        )
        groups.setdefault(sig, ArchitectureGroup(signature=sig)).members.append(pid)

    background = source.random_loci(config.background_loci, seed=7,
                                    window_bp=config.neighbourhood_bp)
    for loc in background:
        calls = labeller.label_genes(
            loc.genes,
            min_z=config.label_min_z,
            min_coverage=config.label_min_coverage,
            cache=label_cache,
        )
        for g in loc.genes:
            call = calls.get(g.gene_id)
            g.label = call.label if call else "unk"
    background_label_sets = [
        {g.label for g in predicted_operon(loc, max_gap=config.operon_max_gap)}
        | ({"__array__"} if loc.arrays else set())
        for loc in background
    ]

    # Partner families are identified by sequence, not by the label "unk":
    # every unlabelled neighbour is called unk, so grouping on the label alone
    # would merge unrelated proteins into one apparent partner family.
    partner_family = _cluster_unknown_partners(loci, novel, config)
    partner_ref = partner_reference_index(source, config.family)

    for sig, grp in groups.items():
        member_loci = [loci[p] for p in grp.members]
        grp.conservation = conservation_across_loci(member_loci)

        records = {p: loci[p].anchor.protein for p in grp.members}
        lineages = {p: loci[p].taxon for p in grp.members}
        spread = taxonomic_spread(records, lineages, identity=config.cluster_identity)
        grp.n_clusters = spread.n_clusters
        grp.n_phyla = spread.n_phyla
        grp.clade_restricted = spread.is_clade_restricted
        grp.spread_summary = spread.summary()

        candidate_label_sets = [
            {g.label for g in predicted_operon(loc, max_gap=config.operon_max_gap)}
            | ({"__array__"} if loc.arrays else set())
            for loc in member_loci
        ]
        tested = sorted(
            {lbl for s in candidate_label_sets for lbl in s} - {config.family}
        )
        raw = {
            lbl: neighbour_enrichment(lbl, candidate_label_sets, background_label_sets)
            for lbl in tested
        }
        # Many neighbour labels are tested per group; without correction some
        # will look enriched by construction.
        flags = benjamini_hochberg([r.p_value for r in raw.values()], alpha=0.05)
        for (lbl, r), sig_flag in zip(raw.items(), flags):
            grp.enrichment[lbl] = {
                "odds_ratio": r.odds_ratio,
                "p_value": r.p_value,
                "significant_bh": bool(sig_flag),
                "summary": r.summary(),
            }
        grp.partner_labels = [
            lbl for lbl in tested
            if grp.conservation.get(lbl, 0.0) >= 0.4 and lbl not in ("__array__",)
        ]

    # -- assemble candidates -------------------------------------------------
    candidates: list[Candidate] = []
    for sig, grp in groups.items():
        for pid in grp.members:
            cand = _build_candidate(
                pid, loci[pid], hit_by_id[pid], assessments[pid], grp,
                partner_family.get(pid, ""), config, profile, partner_ref,
            )
            candidates.append(cand)

    candidates.sort(key=lambda c: (-c.novelty, -float(c.metric("partner_conservation", 0.0))))
    if len(candidates) > config.max_candidates:
        rejected["below the candidate cap after ranking"] += len(candidates) - config.max_candidates
        candidates = candidates[: config.max_candidates]

    result = SearchResult(
        scanned=len(index),
        funnel=funnel,
        profile_hits=len(hits),
        loci=len(loci),
        novel_loci=len(novel),
        groups=[asdict(g) for g in groups.values()],
        candidates=len(candidates),
        rejected=dict(rejected),
        notes=[
            funnel,
            f"{len(loci)} loci assembled; {len(novel)} scored novelty >= {config.min_novelty}",
            f"{len(groups)} distinct architectures among the novel loci",
            f"{len(background)} background loci sampled for enrichment tests",
        ],
    )

    run.write_json("search.json", asdict(result))
    run.write_json("candidates.json", [c.model_dump() for c in candidates])
    run.record_stage(StageRecord(
        stage="search",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"min_z": config.profile_min_z, "min_novelty": config.min_novelty},
        outputs={"hits": len(hits), "loci": len(loci), "novel": len(novel),
                 "groups": len(groups), "candidates": len(candidates)},
        notes=result.notes,
    ))
    return result, candidates


def _cluster_unknown_partners(
    loci: dict[str, Locus],
    novel: dict,
    config: CampaignConfig,
) -> dict[str, str]:
    """Group unlabelled neighbours into families by sequence similarity.

    Two loci both carrying an "unknown partner" only support each other if the
    two unknowns are the same protein. Labelling both ``unk`` and counting them
    together is the mistake this prevents.
    """
    from ..bio.cluster import greedy_cluster

    reps: dict[str, str] = {}
    for pid in novel:
        loc = loci[pid]
        unknowns = [
            g for g in predicted_operon(loc, max_gap=config.operon_max_gap)
            if g.label == "unk" and g.gene_id != loc.anchor.gene_id and g.protein
        ]
        if unknowns:
            reps[pid] = max(unknowns, key=lambda g: len(g.protein)).protein

    if not reps:
        return {}
    clusters = greedy_cluster(reps, identity=0.35, coverage=0.5, prefilter_jaccard=0.15)
    out: dict[str, str] = {}
    for i, c in enumerate(clusters):
        name = f"unkfam{i + 1:02d}"
        for pid in c.members:
            out[pid] = name
    return out


def _build_candidate(
    pid: str,
    locus: Locus,
    hit: Hit,
    assessment,
    group: ArchitectureGroup,
    partner_family: str,
    config: CampaignConfig,
    profile: Profile,
    partner_ref_index: SequenceIndex | None,
) -> Candidate:
    """Turn one novel locus into a candidate with its measurements attached."""
    anchor = locus.anchor
    triad_ok, triad_why = catalytic_triad_intact(anchor.protein)
    low_c = low_complexity_fraction(anchor.protein)
    arr_q, arr_why, arr = _array_quality(locus)
    arr_gap = _operon_distance(locus, arr)

    operon = predicted_operon(locus, max_gap=config.operon_max_gap)
    partners = [g for g in operon if g.gene_id != anchor.gene_id]
    partner_label = partner_family or next(
        (g.label for g in partners if g.label != "unk"), "unk"
    )
    conservation = max(
        [group.conservation.get(g.label, 0.0) for g in partners] or [0.0]
    )
    if partner_family:
        # Conservation of the *family*, not of the "unk" label: how many loci in
        # this architecture group carry a partner from the same cluster.
        conservation = max(conservation, group.conservation.get("unk", 0.0))

    enr = group.enrichment.get("unk") or group.enrichment.get(partner_label) or {}
    p_value = float(enr.get("p_value", 1.0))
    import math
    neglog_p = -math.log10(max(p_value, 1e-300))

    # The sensitive re-check that catches a diverged member of a described
    # family posing as an unknown partner. The strict labelling threshold
    # deliberately misses these; the gate in the rubric reads this number.
    known_e, known_bits = _partner_known_homology(partners, partner_ref_index)

    metrics = {
        "triad_intact": triad_ok,
        "core_coverage": round(hit.core_coverage, 4),
        "aligned_fraction": round(hit.alignment.aligned_len / max(1, profile.length), 4),
        "low_complexity": round(low_c, 4),
        "evalue": hit.evalue,
        "profile_z": round(hit.z, 3),
        "bits": round(hit.bits, 2),
        "novelty": round(assessment.novelty, 4),
        "partner_conservation": round(conservation, 4),
        "partner_known_neglog_e": round(known_e, 2),
        "partner_known_homology_bits": round(known_bits, 1),
        "n_independent_loci": len(group.members),
        "n_clusters": group.n_clusters,
        "n_phyla": group.n_phyla,
        "clade_restricted": group.clade_restricted,
        "association_neglog_p": round(neglog_p, 3),
        "association_significant": bool(enr.get("significant_bh", False)),
        "array_quality": round(arr_q, 4),
        "array_copies": arr.copy_number if arr else 0,
        "array_repeat_len": arr.repeat_length if arr else 0,
        "array_spacer_identity": round(arr.spacer_identity, 3) if arr else 0.0,
        "array_operon_gap_bp": arr_gap,
        "anchor_aa": len(anchor.protein),
        "n_predictions": 0,   # filled in by the report stage
    }

    def sfor(key: str, fallback: Strength = Strength.SUGGESTIVE) -> Strength:
        """Strength from the shared rules, so reports and reviews agree."""
        return strength_for(key, metrics) or fallback

    evidence = [
        Evidence(
            key="catalytic_core", statement=triad_why, value=triad_ok,
            strength=sfor("catalytic_core"),
            method="PROSITE-style motif model of the RT palm domain (motifs A and C)",
            records=[anchor.gene_id], supports=triad_ok,
        ),
        Evidence(
            key="family_membership",
            statement=(
                f"profile hit at {hit.bits:.0f} bits (E={hit.evalue:.1e}, "
                f"z={hit.z:.0f} above the shuffled-decoy null), covering "
                f"{hit.core_coverage:.0%} of the profile's core columns"
            ),
            value=hit.bits, strength=sfor("family_membership"),
            method=f"calibrated PSSM of {profile.length} columns vs. the protein database",
            records=[anchor.gene_id],
        ),
        Evidence(
            key="architecture",
            statement=(
                f"locus architecture is {core_architecture(locus, max_gap=config.operon_max_gap)}; "
                f"{assessment.rationale}"
            ),
            value=assessment.novelty, strength=sfor("architecture"),
            method="operon prediction plus architectural comparison to the described-system catalogue",
            records=[locus.contig_id],
        ),
        Evidence(
            key="partner_association",
            statement=(
                f"the partner travels with the enzyme at {conservation:.0%} of the "
                f"{len(group.members)} loci sharing this architecture; "
                f"{enr.get('summary', 'no enrichment test was possible')}"
                + ("" if enr.get("significant_bh") else
                   " (not significant after Benjamini-Hochberg correction)")
            ),
            value=conservation,
            strength=sfor("partner_association"),
            method="Fisher exact test against background loci from the same genomes, BH-corrected",
            records=[g.gene_id for g in partners][:5],
            supports=conservation > 0.0,
        ),
        Evidence(
            key="noncoding_component", statement=arr_why, value=arr_q,
            strength=sfor("noncoding_component"),
            method="repeat-spacer array detection on the contig, with spacer-distinctness scoring",
            records=[locus.contig_id], supports=arr_q > 0.0,
        ),
        Evidence(
            key="taxonomic_spread", statement=group.spread_summary,
            value=group.n_phyla,
            strength=sfor("taxonomic_spread"),
            method="greedy clustering at 90% identity plus lineage tabulation",
            records=group.members[:5],
            supports=not group.clade_restricted,
        ),
    ]

    if arr is not None and arr_gap > 3000:
        evidence.append(Evidence(
            key="array_placement",
            statement=(
                f"the nearest array is {arr_gap} bp from the predicted operon, "
                f"too far to be part of the same transcriptional unit"
            ),
            value=arr_gap, strength=Strength.STRONG,
            method="distance from the predicted operon bounds to the array",
            records=[locus.contig_id], supports=False,
        ))

    if known_e >= config.partner_known_neglog_e_cut:
        evidence.append(Evidence(
            key="partner_is_known",
            statement=(
                f"the unlabelled partner hits a described family at "
                f"E=1e-{known_e:.0f} ({known_bits:.0f} bits) under a sensitive "
                f"search: it is a diverged member of that family, not an unknown"
            ),
            value=known_e, strength=Strength.STRONG,
            method="sensitive reciprocal search of the partner against the reference families",
            records=[g.gene_id for g in partners][:3], supports=False,
        ))

    return Candidate(
        candidate_id=f"cand-{pid}",
        anchor_protein=pid,
        contig_id=locus.contig_id,
        lineage=locus.taxon,
        architecture=architecture(locus),
        core_architecture=core_architecture(locus, max_gap=config.operon_max_gap),
        group_signature=group.signature,
        novelty=assessment.novelty,
        closest_system=assessment.best_match.name if assessment.best_match else "",
        closest_system_score=assessment.best_match.score if assessment.best_match else 0.0,
        n_family_members=len(group.members),
        metrics=metrics,
        evidence=evidence,
        homolog_loci=group.members[:20],
    )


def _partner_known_homology(
    partners, ref_index: SequenceIndex | None
) -> tuple[float, float]:
    """Best sensitive hit of any partner to the described reference families.

    Returns ``(-log10(E), bits)`` for the strongest hit, or ``(0, 0)`` when
    there is none. The E-value is what the rubric gates on: it is comparable
    across partner families of different lengths, which a bit score is not.
    """
    import math

    if ref_index is None or not partners:
        return 0.0, 0.0
    best_e, best_bits = float("inf"), 0.0
    for g in partners:
        if not g.protein:
            continue
        hits = ref_index.search(g.protein, max_evalue=1.0, limit=1)
        if hits and hits[0].evalue < best_e:
            best_e, best_bits = hits[0].evalue, hits[0].bits
    if best_e == float("inf") or best_e <= 0:
        return (0.0 if best_e == float("inf") else 300.0), best_bits
    return max(0.0, -math.log10(best_e)), best_bits


def partner_reference_index(source: DataSource, family: str) -> SequenceIndex:
    """Reference index for the sensitive partner re-check.

    The anchor family itself is excluded: an RT beside an RT is not evidence
    that the partner is a known protein.
    """
    idx = SequenceIndex()
    for fam, seq in source.reference_proteins().items():
        if fam == family:
            continue
        idx.add(f"ref:{fam}", seq)
    return idx
