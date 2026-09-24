"""Iterative profile construction, PSI-BLAST style.

A single reference sequence is a weak detector. Below roughly 30% identity it
misses most of a family, and no threshold fixes that: the signal is not there
in one sequence. The standard remedy is to bootstrap a profile -- search with
the sequence, build a profile from the confident hits, search again with the
profile, repeat -- and it is the remedy used here.

This matters for gene labelling specifically. A neighbouring gene that goes
unlabelled is read downstream as "unknown partner", which is the observation
the whole campaign is hunting for. A labeller that is merely insensitive
therefore does not fail quietly; it manufactures novelty. The methods-check
stage caught this pipeline doing exactly that on short partner families, which
is what this module is here to fix.

Nothing here consults family annotations. A family model is built from one
reference sequence plus whatever the database itself yields, so the procedure
works identically on a real database where no annotation exists.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .profile import Profile, profile_from_unaligned
from .search import SequenceIndex


@dataclass
class FamilyModel:
    """A bootstrapped model of one described family."""

    family: str
    profile: Profile
    n_members: int
    rounds: int
    member_ids: list[str] = field(default_factory=list)
    converged: bool = False

    def score(self, seq: str, *, search_space: tuple[int, int] | None = None):
        return self.profile.search(seq, search_space=search_space)


def bootstrap_family_model(
    family: str,
    reference: str,
    index: SequenceIndex,
    *,
    rounds: int = 2,
    inclusion_evalue: float = 1e-6,
    max_members: int = 30,
    min_members: int = 3,
    calibration_n: int = 50,
) -> FamilyModel:
    """Build a family profile from one reference sequence and a database.

    Each round adds the hits that clear ``inclusion_evalue`` and rebuilds the
    profile from them. Inclusion is deliberately strict: profile bootstrapping
    drifts when a marginal hit is admitted early, and a model that has drifted
    into a neighbouring family will label genes confidently and wrongly.
    """
    # Round 1 is a plain sequence search, so no profile is built before there
    # is anything to build it from. Each subsequent round rebuilds from the
    # hits the previous model found. Profile construction includes decoy
    # calibration and is the expensive step, so it happens once per round and
    # never speculatively.
    hits = index.search(reference, max_evalue=inclusion_evalue, limit=max_members)
    member_ids = [h.target_id for h in hits]
    members = [reference] + [index.get(i) for i in member_ids]
    profile = profile_from_unaligned(
        family, members[:max_members], calibration_n=calibration_n
    )
    converged = False
    used_rounds = 1

    for _ in range(rounds - 1):
        hits = index.profile_search(
            profile,
            max_evalue=inclusion_evalue,
            min_z=8.0,
            min_core_coverage=0.5,
            limit=max_members,
        )
        new_ids = [h.target_id for h in hits]
        if set(new_ids) == set(member_ids):
            converged = True
            break
        used_rounds += 1
        member_ids = new_ids
        members = [reference] + [index.get(i) for i in new_ids]
        if len(members) < min_members:
            # Too few confident relatives to build a profile worth having;
            # keep the previous model rather than fitting noise.
            break
        profile = profile_from_unaligned(
            family, members[:max_members], calibration_n=calibration_n
        )

    return FamilyModel(
        family=family,
        profile=profile,
        n_members=len(members),
        rounds=used_rounds,
        member_ids=member_ids,
        converged=converged,
    )


def bootstrap_all(
    references: dict[str, str],
    index: SequenceIndex,
    *,
    rounds: int = 2,
    inclusion_evalue: float = 1e-6,
    max_members: int = 30,
    calibration_n: int = 50,
    progress=None,
) -> dict[str, FamilyModel]:
    """Bootstrap a model for every reference family."""
    models: dict[str, FamilyModel] = {}
    for i, (family, seq) in enumerate(sorted(references.items())):
        if progress is not None:
            progress(i, len(references), family)
        models[family] = bootstrap_family_model(
            family, seq, index,
            rounds=rounds, inclusion_evalue=inclusion_evalue,
            max_members=max_members, calibration_n=calibration_n,
        )
    return models
