"""Profile models: build a PSSM from a seed alignment, calibrate it, scan with it.

This is the studio's stand-in for ``hmmbuild``/``hmmsearch``. It is a
position-specific scoring matrix with sequence weighting and Dirichlet-style
pseudocounts rather than a full profile HMM -- no per-position insert/delete
states -- but it carries the two properties the workflow depends on: scores are
log-odds against an explicit background, and the score-to-significance mapping
is calibrated against a shuffled decoy set rather than assumed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .align import Alignment, encode_batch, profile_align, profile_scores_batch
from .seq import AA, AA_INDEX, BACKGROUND_AA, clean_protein

_BG = np.array([BACKGROUND_AA[a] for a in AA], dtype=np.float64)


@dataclass
class Profile:
    """A calibrated position-specific scoring matrix."""

    name: str
    pssm: np.ndarray                      # (L, 21) log-odds in bits
    consensus: str
    information: np.ndarray               # (L,) bits of information per column
    n_seed: int
    # Calibration against reversed-sequence decoys: score -> significance.
    decoy_mu: float = 0.0
    decoy_sigma: float = 1.0
    core_columns: list[int] = field(default_factory=list)

    @property
    def length(self) -> int:
        return int(self.pssm.shape[0])

    def z(self, score: float) -> float:
        """Standard deviations above the decoy mean."""
        return (score - self.decoy_mu) / (self.decoy_sigma or 1.0)

    def search(self, target: str, *, search_space: tuple[int, int] | None = None) -> Alignment:
        return profile_align(self.pssm, clean_protein(target), search_space=search_space)


def _sequence_weights(cols: list[list[str]], n: int) -> np.ndarray:
    """Henikoff position-based sequence weights.

    Down-weights groups of near-identical seed sequences so that a profile built
    from a lopsided seed is not dominated by whichever clade happened to be
    sequenced most.
    """
    w = np.zeros(n, dtype=np.float64)
    for col in cols:
        present = [(k, c) for k, c in enumerate(col) if c != "-"]
        if not present:
            continue
        counts: dict[str, int] = {}
        for _, c in present:
            counts[c] = counts.get(c, 0) + 1
        r = len(counts)
        for k, c in present:
            w[k] += 1.0 / (r * counts[c])
    if w.sum() == 0:
        return np.ones(n, dtype=np.float64)
    # Scale to sum to the sequence count, not to 1: the weights stand in for
    # observation counts, so a fixed pseudocount mass must be small relative to
    # a deep column and comparable to a shallow one.
    return w * (n / w.sum())


def build_profile(
    name: str,
    seed: list[str],
    *,
    max_gap_fraction: float = 0.5,
    pseudocount: float = 1.5,
    core_information: float = 1.5,
    calibration_n: int = 200,
    rng: np.random.Generator | None = None,
) -> Profile:
    """Build a PSSM from a seed alignment (equal-length, ``-`` for gaps).

    Columns with more than ``max_gap_fraction`` gaps are dropped, which turns a
    padded alignment into a match-state-only model. Columns whose information
    content exceeds ``core_information`` bits are recorded as ``core_columns``:
    the search stage requires a hit to cover them before it will call the domain
    present, so a fragment that clips only the variable edges of a domain does
    not become a candidate.
    """
    if not seed:
        raise ValueError("seed alignment is empty")
    width = len(seed[0])
    if any(len(s) != width for s in seed):
        raise ValueError("seed alignment rows must be equal length")

    rows = [s.upper() for s in seed]
    cols = [[r[j] for r in rows] for j in range(width)]
    keep = [j for j, col in enumerate(cols)
            if sum(1 for c in col if c == "-") / len(col) <= max_gap_fraction]
    if not keep:
        raise ValueError("no columns survived the gap filter")

    weights = _sequence_weights([cols[j] for j in keep], len(rows))

    pssm = np.zeros((len(keep), 21), dtype=np.float32)
    info = np.zeros(len(keep), dtype=np.float64)
    consensus: list[str] = []

    for out_j, j in enumerate(keep):
        obs = np.zeros(20, dtype=np.float64)
        for k, c in enumerate(cols[j]):
            idx = AA_INDEX.get(c)
            if idx is not None and idx < 20:
                obs[idx] += weights[k]
        total = obs.sum()
        # Mix observed counts with background pseudocounts; the pseudocount mass
        # is fixed rather than proportional so shallow columns stay conservative.
        freq = (obs + pseudocount * _BG) / (total + pseudocount)
        odds = np.log2(freq / _BG)
        pssm[out_j, :20] = odds
        pssm[out_j, 20] = -1.0  # X: mild penalty, never informative
        info[out_j] = float(np.sum(freq * np.log2(freq / _BG)))
        consensus.append(AA[int(np.argmax(freq))])

    core = [j for j in range(len(keep)) if info[j] >= core_information]

    prof = Profile(
        name=name,
        pssm=pssm,
        consensus="".join(consensus),
        information=info,
        n_seed=len(rows),
        core_columns=core,
    )
    _calibrate(prof, rows, calibration_n, rng or np.random.default_rng(0))
    return prof


def _calibrate(prof: Profile, rows: list[str], n: int, rng: np.random.Generator) -> None:
    """Estimate the null score distribution from shuffled decoys.

    Shuffling preserves each seed sequence's amino-acid composition while
    destroying its order, so the resulting mean and spread measure how much
    score this profile hands out for composition alone. Every reported hit is
    quoted as a z-score against this null.
    """
    pool = [clean_protein(r.replace("-", "")) for r in rows if r.strip("-")]
    if not pool:
        prof.decoy_mu, prof.decoy_sigma = 0.0, 1.0
        return
    decoys: list[str] = []
    for i in range(n):
        src = list(pool[i % len(pool)])
        rng.shuffle(src)
        decoys.append("".join(src))
    # Batched: calibration is n alignments per profile and a campaign builds
    # one profile per family, so this is the difference between seconds and
    # minutes of start-up.
    enc, lengths = encode_batch(decoys)
    arr = profile_scores_batch(prof.pssm, enc, lengths).astype(np.float64)
    prof.decoy_mu = float(arr.mean())
    prof.decoy_sigma = float(arr.std() or 1.0)


def profile_from_unaligned(
    name: str,
    seqs: list[str],
    **kwargs,
) -> Profile:
    """Build a profile from unaligned sequences via progressive anchor alignment.

    The longest sequence becomes the anchor; every other sequence is aligned to
    it semi-globally and projected onto anchor coordinates. This is cruder than
    a true progressive MSA but adequate for seeds drawn from one domain family,
    and it keeps the studio free of external alignment binaries.
    """
    from .align import align

    cleaned = [clean_protein(s) for s in seqs if s.strip()]
    if not cleaned:
        raise ValueError("no sequences to align")
    anchor = max(cleaned, key=len)
    rows = [anchor]
    for s in cleaned:
        if s is anchor:
            continue
        a = align(anchor, s, local=False)
        row = ["-"] * len(anchor)
        qi = a.q_start
        for qc, tc in zip(a.q_aln, a.t_aln):
            if qc != "-":
                if tc != "-" and qi < len(anchor):
                    row[qi] = tc
                qi += 1
        rows.append("".join(row))
    return build_profile(name, rows, **kwargs)
