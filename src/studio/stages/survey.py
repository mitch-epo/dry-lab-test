"""Stage 1 -- survey: read the literature, build the family model.

"One typical pattern begins with a survey of a given protein family. Claude
reads the relevant literature..."

Two things come out of this stage and everything downstream depends on both:

* a **family profile**, calibrated against a decoy null so that later hits can
  be quoted as a z-score rather than as an unanchored number; and
* a **catalogue of described systems**, turned from prose into the machine-
  checkable architectural predicates the search stage measures novelty
  against.

The second is the part that is easy to skip and expensive to skip. Without an
explicit, falsifiable statement of what is already known, "fits no described
system" has no meaning, and a survey will confidently rediscover retrons.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from typing import Any

from ..bio.profile import Profile, profile_from_unaligned
from ..bio.search import SequenceIndex
from ..campaign import CampaignConfig, Run, load_catalogue
from ..schemas import StageRecord
from pathlib import Path

from ..sources.base import DataSource, Publication


@dataclass
class SurveyResult:
    """What the survey stage establishes."""

    family: str
    profile_length: int
    profile_core_columns: int
    profile_seed_n: int
    decoy_mu: float
    decoy_sigma: float
    catalogue_size: int
    described_systems: list[dict] = field(default_factory=list)
    publications: list[dict] = field(default_factory=list)
    discriminating_features: list[str] = field(default_factory=list)
    known_pitfalls: list[str] = field(default_factory=list)
    reference_families: list[str] = field(default_factory=list)
    database_size: int = 0
    family_models: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def _discriminating_features(systems) -> list[str]:
    """Turn the catalogue into the specific observations that separate systems.

    This is the survey's most useful product: a written statement of what would
    have to be true for a locus to be each described system, and therefore what
    observation would rule each one out.
    """
    out: list[str] = []
    for s in systems:
        req = " + ".join(s.required_labels) or "(none)"
        bits = [f"{s.name}: requires {req}"]
        if s.requires_array:
            rng = s.array_repeat_range
            bits.append(
                f"requires a repeat array{f' of {rng[0]}-{rng[1]} bp repeats' if rng else ''}"
            )
        if s.forbids_array:
            bits.append("has NO multi-copy repeat array (an array rules it out)")
        if s.forbidden_labels:
            bits.append(f"never carries {', '.join(s.forbidden_labels)}")
        if s.typical_anchor_aa:
            bits.append(f"anchor typically {s.typical_anchor_aa[0]}-{s.typical_anchor_aa[1]} aa")
        out.append("; ".join(bits))
    return out


def run_survey(
    run: Run,
    config: CampaignConfig,
    source: DataSource,
    *,
    seed_limit: int = 14,
) -> tuple[SurveyResult, Profile, SequenceIndex, list, "FamilyLabeller"]:
    """Execute the survey. Returns the result plus the objects later stages need."""
    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    systems, _targets, doc = load_catalogue(config.catalogue, base=run.dir.parent.parent)
    pubs: list[Publication] = source.literature(
        f"{config.family} reverse transcriptase retron CRISPR "
        f"diversity generating retroelement gene neighbourhood",
        limit=20,
    )

    seed = source.family_seed(config.family, limit=seed_limit)
    if len(seed) < 3:
        raise RuntimeError(
            f"only {len(seed)} seed sequences for family {config.family!r}: "
            "a profile built from this would not be calibratable"
        )
    profile = profile_from_unaligned(config.family, seed)

    index = SequenceIndex()
    for rec in source.iter_proteins():
        index.add(rec.protein_id, rec.sequence)

    pitfalls: list[str] = []
    for p in pubs:
        if "pitfall" in p.pub_id or "failure" in p.title.lower():
            pitfalls.extend(p.claims)

    labeller, cached = load_or_build_labeller(
        source, index, run.dir.parent.parent / ".cache"
    )

    result = SurveyResult(
        family=config.family,
        profile_length=profile.length,
        profile_core_columns=len(profile.core_columns),
        profile_seed_n=len(seed),
        decoy_mu=round(profile.decoy_mu, 3),
        decoy_sigma=round(profile.decoy_sigma, 3),
        catalogue_size=len(systems),
        described_systems=[asdict(s) for s in systems],
        publications=[
            {"id": p.pub_id, "title": p.title, "claims": p.claims, "source": p.source}
            for p in pubs
        ],
        discriminating_features=_discriminating_features(systems),
        known_pitfalls=pitfalls,
        reference_families=sorted(source.reference_proteins()),
        database_size=len(index),
        notes=[
            f"catalogue v{doc.get('catalogue_version', '?')} with "
            f"{len(systems)} described systems",
            f"profile built from {len(seed)} seed sequences, "
            f"{len(profile.core_columns)} core columns above 1.5 bits",
            f"decoy calibration: mu={profile.decoy_mu:.1f} sigma={profile.decoy_sigma:.1f} "
            f"from shuffled seeds",
            f"{'loaded' if cached else 'bootstrapped'} {len(labeller.models)} "
            f"family profiles for gene labelling",
        ],
    )
    result.family_models = labeller.summary()

    run.write_json("survey.json", asdict(result))
    run.record_stage(StageRecord(
        stage="survey",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"family": config.family, "catalogue": config.catalogue},
        outputs={
            "profile_length": profile.length,
            "database_size": len(index),
            "catalogue_size": len(systems),
            "publications": len(pubs),
            "family_models": len(labeller.models),
        },
        notes=result.notes,
    ))
    return result, profile, index, systems, labeller


@dataclass
class LabelCall:
    """The outcome of labelling one gene."""

    label: str
    z: float = 0.0
    evalue: float = float("inf")
    bits: float = 0.0
    coverage: float = 0.0
    runner_up: str = ""
    runner_up_z: float = 0.0

    @property
    def neglog_e(self) -> float:
        """-log10(E), with 0 meaning no hit. Convenient in rubric expressions."""
        import math

        if self.evalue == float("inf"):
            return 0.0
        if self.evalue <= 0:
            return 300.0
        return max(0.0, -math.log10(self.evalue))


@dataclass
class _Candidate:
    """One family's scoring of one gene, before any threshold is applied.

    ``coverage`` is -1 until a traceback has been run for this pair. Scanning
    stores only the score, which needs no traceback; coverage is resolved on
    demand for the handful of pairs a caller actually asks about.
    """

    family: str
    z: float
    raw: float = 0.0
    evalue: float = float("inf")
    bits: float = 0.0
    coverage: float = -1.0


class FamilyLabeller:
    """Labels genes by scoring them against bootstrapped family profiles.

    WHY PROFILES. A single reference sequence cannot detect its own family
    below roughly 30% identity, so a sequence-based labeller leaves diverged
    members of *described* families unlabelled. Downstream, an unlabelled
    neighbour reads as "unknown partner" -- precisely the observation this
    campaign treats as novel. An insensitive labeller therefore does not lose
    candidates, it invents them. The methods-check stage caught this pipeline
    doing exactly that to short partner families.

    WHY Z-SCORES. Significance is measured against each profile's own
    shuffled-decoy null, which is comparable across families of very different
    lengths. Neither a raw bit score nor a single fixed E-value cut-off is.

    WHY SCAN ONCE. The obvious implementation -- score each gene against every
    family profile on demand -- costs one full alignment per (gene, family)
    pair and is hopeless at campaign scale. Instead each family model is run
    over the whole database once through the k-mer and ungapped funnel, and
    labelling afterwards is a dictionary lookup. Both threshold tiers read the
    same scan, so the strict and sensitive passes cost nothing extra.
    """

    def __init__(
        self,
        models: dict,
        scores: dict[str, list[_Candidate]] | None = None,
        *,
        index: SequenceIndex | None = None,
        database_residues: int = 0,
    ):
        self.models = models
        self.scores = scores or {}
        self.index = index
        self.database_residues = database_residues

    def _resolve(self, gene_id: str, cand: _Candidate) -> _Candidate:
        """Fill in coverage and E-value for one pair, running the traceback once."""
        if cand.coverage >= 0.0 or self.index is None:
            return cand
        model = self.models[cand.family]
        seq = self.index.get(gene_id)
        a = model.profile.search(
            seq, search_space=(len(seq), self.database_residues or len(seq))
        )
        cand.evalue = a.evalue
        cand.bits = a.bits
        cand.coverage = a.aligned_len / max(1, min(model.profile.length, len(seq)))
        return cand

    @classmethod
    def build(
        cls,
        source: DataSource,
        index: SequenceIndex,
        *,
        scan_min_z: float = 5.0,
        progress=None,
    ) -> "FamilyLabeller":
        """Bootstrap every family model, then scan each over the database.

        ``scan_min_z`` is the floor at which a (gene, family) score is
        retained. It sits below every threshold a caller will use, so the
        stored scan can answer both the strict and the sensitive query.
        """
        from ..bio.iterate import bootstrap_all

        refs = source.reference_proteins()
        models = bootstrap_all(refs, index, progress=progress)

        scores: dict[str, list[_Candidate]] = {}
        for i, (family, model) in enumerate(sorted(models.items())):
            if progress is not None:
                progress(len(models) + i, 2 * len(models), f"scan:{family}")
            for target_id, z, raw in index.profile_scan(
                model.profile, min_z=scan_min_z
            ):
                scores.setdefault(target_id, []).append(_Candidate(family, z, raw))
        for gid in scores:
            scores[gid].sort(key=lambda c: -c.z)
        return cls(models, scores, index=index, database_residues=index.total_residues)

    def summary(self) -> list[str]:
        return [
            f"{m.family}: profile of {m.profile.length} columns from "
            f"{m.n_members} members after {m.rounds} round(s)"
            f"{' (converged)' if m.converged else ''}"
            for m in sorted(self.models.values(), key=lambda m: m.family)
        ]

    def call_by_id(self, gene_id: str, *, min_z: float, min_coverage: float) -> LabelCall:
        """Best family call for a gene already present in the scanned database."""
        ranked = self.scores.get(gene_id, [])
        if not ranked:
            return LabelCall("unk")
        best = self._resolve(gene_id, ranked[0])
        runner = ranked[1] if len(ranked) > 1 else None
        if best.z >= min_z and best.coverage >= min_coverage:
            return LabelCall(
                best.family, z=best.z, evalue=best.evalue, bits=best.bits,
                coverage=best.coverage,
                runner_up=runner.family if runner else "",
                runner_up_z=runner.z if runner else 0.0,
            )
        return LabelCall(
            "unk", z=best.z, evalue=best.evalue, bits=best.bits,
            coverage=best.coverage, runner_up=best.family, runner_up_z=best.z,
        )

    def call(self, protein: str, *, min_z: float, min_coverage: float) -> LabelCall:
        """Score an arbitrary protein against every family model.

        The slow path, for sequences that are not in the scanned database.
        """
        if not protein:
            return LabelCall("unk")
        space = (len(protein), self.database_residues or len(protein))
        ranked: list[_Candidate] = []
        for family, model in self.models.items():
            a = model.profile.search(protein, search_space=space)
            ranked.append(_Candidate(
                family, model.profile.z(a.score), a.score, a.evalue, a.bits,
                a.aligned_len / max(1, min(model.profile.length, len(protein))),
            ))
        ranked.sort(key=lambda c: -c.z)
        best, runner = ranked[0], (ranked[1] if len(ranked) > 1 else None)
        if best.z >= min_z and best.coverage >= min_coverage:
            return LabelCall(
                best.family, z=best.z, evalue=best.evalue, bits=best.bits,
                coverage=best.coverage,
                runner_up=runner.family if runner else "",
                runner_up_z=runner.z if runner else 0.0,
            )
        return LabelCall(
            "unk", z=best.z, evalue=best.evalue, bits=best.bits,
            coverage=best.coverage, runner_up=best.family, runner_up_z=best.z,
        )

    def label_genes(
        self,
        genes,
        *,
        min_z: float,
        min_coverage: float,
        cache: dict[str, LabelCall] | None = None,
    ) -> dict[str, LabelCall]:
        """Label every gene from the stored scan, falling back to a full score."""
        out: dict[str, LabelCall] = {}
        for g in genes:
            if cache is not None and g.gene_id in cache:
                out[g.gene_id] = cache[g.gene_id]
                continue
            if g.gene_id in self.scores:
                call = self.call_by_id(g.gene_id, min_z=min_z, min_coverage=min_coverage)
            else:
                call = self.call(g.protein, min_z=min_z, min_coverage=min_coverage)
            out[g.gene_id] = call
            if cache is not None:
                cache[g.gene_id] = call
        return out


def build_reference_index(source: DataSource) -> tuple[SequenceIndex, dict[str, str]]:
    """Index the reference proteins themselves, for single-sequence lookups.

    Used where a plain sequence search is what is wanted -- the sensitive
    partner re-check in triage, for instance. Gene labelling goes through
    :class:`FamilyLabeller` instead, for the reasons given there.
    """
    refs = source.reference_proteins()
    idx = SequenceIndex()
    label_of: dict[str, str] = {}
    for family, seq in refs.items():
        rid = f"ref:{family}"
        idx.add(rid, seq)
        label_of[rid] = family
    return idx, label_of


# --- labeller persistence ----------------------------------------------------
#
# Bootstrapping two dozen family profiles is minutes of work and depends only
# on the reference set and the database, neither of which changes between runs
# of the same campaign. Caching it keeps the iteration loop short; the key
# covers both inputs so a changed corpus invalidates the cache rather than
# silently reusing models built for different data.

def _labeller_key(source: DataSource, index: SequenceIndex) -> str:
    import hashlib

    h = hashlib.sha256()
    for family, seq in sorted(source.reference_proteins().items()):
        h.update(family.encode())
        h.update(seq.encode())
    h.update(str(len(index)).encode())
    for sid in index.ids:
        h.update(sid.encode())
    return h.hexdigest()[:20]


def load_or_build_labeller(
    source: DataSource,
    index: SequenceIndex,
    cache_dir: Path | None,
    *,
    progress=None,
) -> tuple["FamilyLabeller", bool]:
    """Return ``(labeller, from_cache)``."""
    import pickle

    if cache_dir is None:
        return FamilyLabeller.build(source, index, progress=progress), False

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"labeller-{_labeller_key(source, index)}.pkl"
    if path.exists():
        try:
            with path.open("rb") as fh:
                models, scores = pickle.load(fh)
            return (
                FamilyLabeller(
                    models, scores, index=index,
                    database_residues=index.total_residues,
                ),
                True,
            )
        except Exception:
            # A cache that cannot be read is replaced, not worked around.
            path.unlink(missing_ok=True)

    labeller = FamilyLabeller.build(source, index, progress=progress)
    # Scores are re-derived from the scan; the resolved coverage fields are
    # dropped so a reload recomputes them against the current index.
    fresh = {
        gid: [_Candidate(c.family, c.z, c.raw) for c in cands]
        for gid, cands in labeller.scores.items()
    }
    with path.open("wb") as fh:
        pickle.dump((labeller.models, fresh), fh)
    return labeller, False
