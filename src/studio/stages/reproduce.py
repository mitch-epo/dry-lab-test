"""Stage 2 -- reproduce: check the methods before trusting them.

"Claude reads the relevant literature and reproduces the established results
from public data to check its methods."

This is the stage that makes the rest of the campaign worth reading, and it is
the one that is easiest to leave out. Before the pipeline is permitted to
propose anything new, it must re-find what is already known, using the same
code path it will use for discovery. If it cannot recover retrons, its opinion
that something is *not* a retron is worthless.

Three checks run here:

1. **Positive controls.** For each described system with a reproduction target,
   an exhaustive reference search establishes the loci that system should
   cover, and the production path is scored against it for recall and
   precision. The two paths are deliberately different -- the reference uses
   single-sequence searches at a permissive threshold, the production path uses
   the calibrated profile at the campaign's strict threshold -- so agreement is
   informative rather than tautological.

2. **Negative control on sequence.** Composition-matched shuffles of real
   proteins must not pass the profile threshold. If they do, the threshold is
   measuring composition.

3. **Negative control on the array detector.** Shuffled contig DNA must not
   yield regularly spaced arrays. If it does, every array claim downstream is
   an artefact.

Failing any of these sets ``gate_passed = False`` and the campaign stops.
"""

from __future__ import annotations

import random
import time
from dataclasses import asdict, dataclass, field

from ..bio.context import Locus, match_known_system, predicted_operon
from ..bio.profile import Profile
from ..bio.repeats import find_repeat_arrays
from ..bio.search import SequenceIndex
from ..campaign import CampaignConfig, Run, load_catalogue
from ..schemas import ReproductionResult, StageRecord
from ..sources.base import DataSource
from .survey import FamilyLabeller


@dataclass
class ControlResult:
    """One negative control."""

    name: str
    description: str
    observed: float
    allowed: float
    passed: bool

    def line(self) -> str:
        flag = "PASS" if self.passed else "FAIL"
        return f"[{flag}] {self.name}: observed {self.observed:.3g}, allowed <= {self.allowed:.3g}"


@dataclass
class ReproduceResult:
    controls: list[ControlResult] = field(default_factory=list)
    systems: list[ReproductionResult] = field(default_factory=list)
    gate_passed: bool = False
    notes: list[str] = field(default_factory=list)

    def report(self) -> str:
        lines = ["Methods check", "=" * 13, ""]
        lines += [c.line() for c in self.controls]
        lines.append("")
        lines += [s.line() for s in self.systems]
        lines += ["", f"GATE: {'PASSED' if self.gate_passed else 'FAILED'}"]
        return "\n".join(lines)


def _anchor_loci_by_profile(
    profile: Profile,
    index: SequenceIndex,
    source: DataSource,
    config: CampaignConfig,
) -> dict[str, Locus]:
    """Production path: profile hits -> loci."""
    hits = index.profile_search(
        profile,
        max_evalue=config.profile_max_evalue,
        min_z=config.profile_min_z,
        min_core_coverage=config.profile_min_core_coverage,
    )
    out: dict[str, Locus] = {}
    for h in hits:
        loc = source.locus(h.target_id, window_bp=config.neighbourhood_bp)
        if loc is not None:
            out[h.target_id] = loc
    return out


def _anchor_loci_by_reference(
    source: DataSource,
    index: SequenceIndex,
    family: str,
    *,
    max_evalue: float,
    window_bp: int,
) -> dict[str, Locus]:
    """Reference path: a permissive single-sequence search for the anchor family."""
    refs = source.reference_proteins()
    query = refs.get(family)
    if not query:
        return {}
    hits = index.search(query, max_evalue=max_evalue, limit=100_000)
    out: dict[str, Locus] = {}
    for h in hits:
        loc = source.locus(h.target_id, window_bp=window_bp)
        if loc is not None:
            out[h.target_id] = loc
    return out


def _apply_labels(
    loci: dict[str, Locus],
    labeller: FamilyLabeller,
    *,
    min_z: float,
    min_coverage: float,
    cache: dict | None = None,
) -> None:
    """Label every gene in every locus, in place.

    The cache is shared across the production and reference passes: the two
    overlap heavily, and a gene's label does not depend on which pass asked.
    """
    for loc in loci.values():
        calls = labeller.label_genes(
            loc.genes, min_z=min_z, min_coverage=min_coverage, cache=cache,
        )
        for g in loc.genes:
            call = calls.get(g.gene_id)
            g.label = call.label if call else "unk"


def _matching(loci: dict[str, Locus], system) -> set[str]:
    """Loci whose architecture matches a described system."""
    out: set[str] = set()
    for pid, loc in loci.items():
        if match_known_system(loc, system).is_match:
            out.add(pid)
    return out


def run_reproduce(
    run: Run,
    config: CampaignConfig,
    source: DataSource,
    profile: Profile,
    index: SequenceIndex,
    systems: list,
    labeller: FamilyLabeller,
    *,
    n_shuffles: int = 150,
    n_dna_shuffles: int = 25,
    seed: int = 11,
) -> ReproduceResult:
    """Execute the methods check. Sets ``run.manifest.gate_passed``."""
    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    rng = random.Random(seed)
    result = ReproduceResult()

    # -- negative control 1: composition-matched protein shuffles -----------
    pool = [r.sequence for r in source.iter_proteins()]
    rng.shuffle(pool)
    false_hits = 0
    for seq in pool[:n_shuffles]:
        chars = list(seq)
        rng.shuffle(chars)
        a = profile.search("".join(chars))
        if a.evalue <= config.profile_max_evalue and profile.z(a.score) >= config.profile_min_z:
            false_hits += 1
    rate = false_hits / max(1, min(n_shuffles, len(pool)))
    result.controls.append(ControlResult(
        name="shuffled-protein-false-positive-rate",
        description=(
            "composition-matched shuffles passing the profile threshold; a "
            "non-zero rate means the threshold is scoring composition"
        ),
        observed=rate, allowed=0.01, passed=rate <= 0.01,
    ))

    # -- negative control 2: shuffled DNA must not produce arrays -----------
    contig_ids = sorted({r.contig_id for r in source.iter_proteins() if r.contig_id})
    rng.shuffle(contig_ids)
    spurious = 0
    checked = 0
    for cid in contig_ids[:n_dna_shuffles]:
        dna = source.contig_dna(cid)
        if not dna:
            continue
        checked += 1
        bases = list(dna)
        rng.shuffle(bases)
        arrays = find_repeat_arrays("".join(bases))
        if any(a.copy_number >= 3 and a.spacer_regularity > 0.75 for a in arrays):
            spurious += 1
    arr_rate = spurious / max(1, checked)
    result.controls.append(ControlResult(
        name="shuffled-dna-array-false-positive-rate",
        description=(
            "regularly spaced arrays called in mononucleotide-shuffled contigs; "
            "any appreciable rate invalidates every array claim downstream"
        ),
        observed=arr_rate, allowed=0.04, passed=arr_rate <= 0.04,
    ))

    # -- positive controls ---------------------------------------------------
    _, targets, _ = load_catalogue(config.catalogue, base=run.dir.parent.parent)

    production = _anchor_loci_by_profile(profile, index, source, config)
    reference = _anchor_loci_by_reference(
        source, index, config.family,
        max_evalue=config.sensitive_label_max_evalue,
        window_bp=config.neighbourhood_bp,
    )
    # Separate caches: the two passes use different thresholds, so a label
    # computed under one is not valid under the other.
    _apply_labels(production, labeller,
                  min_z=config.label_min_z,
                  min_coverage=config.label_min_coverage, cache={})
    _apply_labels(reference, labeller,
                  min_z=config.sensitive_label_min_z,
                  min_coverage=config.sensitive_label_min_coverage, cache={})

    by_id = {s.system_id: s for s in systems}
    for target in targets:
        sid = target["system_id"]
        system = by_id.get(sid)
        if system is None:
            continue
        expected = _matching(reference, system)
        recovered = _matching(production, system)
        tp = len(expected & recovered)
        fp = len(recovered - expected)
        result.systems.append(ReproductionResult(
            system_id=sid,
            expected=len(expected),
            recovered=tp,
            false_positives=fp,
            recall=tp / len(expected) if expected else 0.0,
            precision=tp / len(recovered) if recovered else 0.0,
            min_recall=float(target.get("min_recall", 0.7)),
            min_precision=float(target.get("min_precision", 0.8)),
        ))

    result.gate_passed = (
        all(c.passed for c in result.controls)
        and all(s.passed for s in result.systems)
        and bool(result.systems)
    )
    result.notes = [
        f"production path recovered {len(production)} anchor loci at z>="
        f"{config.profile_min_z}, E<={config.profile_max_evalue:g}",
        f"reference path recovered {len(reference)} anchor loci at "
        f"E<={config.sensitive_label_max_evalue:g}",
        f"labels: strict z>={config.label_min_z}, sensitive z>={config.sensitive_label_min_z}",
        f"{len(result.systems)} positive controls, {len(result.controls)} negative controls",
    ]
    if not result.systems:
        result.notes.append(
            "no positive controls were evaluated: the catalogue declares no "
            "reproduction targets, so the gate cannot pass"
        )

    run.write_json("reproduce.json", asdict(result))
    run.write_text("reproduce.txt", result.report() + "\n")
    run.manifest.gate_passed = result.gate_passed
    run.record_stage(StageRecord(
        stage="reproduce",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"targets": [t["system_id"] for t in targets]},
        outputs={
            "gate_passed": result.gate_passed,
            "controls_passed": sum(1 for c in result.controls if c.passed),
            "systems_passed": sum(1 for s in result.systems if s.passed),
        },
        notes=result.notes,
    ))
    return result
