"""Detection of tandem repeat arrays in DNA.

The system at the centre of the article is defined partly by a non-coding
component: an array of evenly spaced DNA repeats sitting beside the enzyme.
Arrays like this are what distinguishes several described RT-associated systems
from one another -- a CRISPR array, a retron's msr/msd, a diversity-generating
retroelement's template repeat -- so the studio treats array architecture as a
first-class, measurable feature rather than as free text.

The detector is CRISPR-array-shaped: find a short unit that recurs at a regular
period, verify the units are more similar to each other than the intervening
spacers are, and report the geometry.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from statistics import mean, pstdev

import numpy as np

from .seq import hamming, revcomp


@dataclass
class RepeatArray:
    """A detected repeat-spacer array."""

    start: int
    end: int
    repeat_consensus: str
    repeat_positions: list[int] = field(default_factory=list)
    spacer_lengths: list[int] = field(default_factory=list)
    spacers: list[str] = field(default_factory=list)
    repeat_identity: float = 0.0     # mean identity of units to the consensus
    spacer_identity: float = 0.0     # mean pairwise identity between spacers
    strand: int = 1

    @property
    def copy_number(self) -> int:
        return len(self.repeat_positions)

    @property
    def repeat_length(self) -> int:
        return len(self.repeat_consensus)

    @property
    def mean_spacer(self) -> float:
        return mean(self.spacer_lengths) if self.spacer_lengths else 0.0

    @property
    def spacer_regularity(self) -> float:
        """1.0 when every spacer is the same length; falls off with the spread.

        Evenly spaced arrays -- the phrase used to describe the system in the
        article -- score near 1. Chance repeats in a genome almost never do.
        """
        if len(self.spacer_lengths) < 2:
            return 0.0
        m = mean(self.spacer_lengths)
        if m == 0:
            return 0.0
        return max(0.0, 1.0 - pstdev(self.spacer_lengths) / m)

    @property
    def length(self) -> int:
        return self.end - self.start

    def classify(self) -> str:
        """A coarse architectural label, not an identification.

        The label is a hypothesis handed to the reporting stage; the triage
        stage is expected to challenge it against the known-system catalogue.
        """
        if self.spacer_identity > 0.6:
            return "tandem-repeat (near-identical spacers; likely structural or VNTR)"
        if 23 <= self.repeat_length <= 50 and 25 <= self.mean_spacer <= 50 \
                and self.copy_number >= 3 and self.spacer_regularity > 0.8:
            return "CRISPR-like (repeat-spacer array with distinct spacers)"
        if self.copy_number >= 3 and self.spacer_regularity > 0.75:
            return "regularly spaced non-coding array"
        return "irregular repeat cluster"

    def class_token(self) -> str:
        """Coarse architectural class, for grouping loci by architecture.

        Copy number and repeat length vary between instances of one real
        system, so a signature carrying them makes every locus unique and
        destroys the recurrence evidence that a shared architecture is
        supposed to provide. The class is the part that is shared; the
        geometry belongs in the report, not in the group key.
        """
        if self.spacer_identity > 0.6:
            return "tandem-repeat"
        if self.copy_number >= 3 and self.spacer_regularity > 0.75:
            return "spacer-array"
        return "irregular-repeat"

    def summary(self) -> str:
        return (
            f"{self.copy_number}x {self.repeat_length} bp repeat, "
            f"spacers {self.mean_spacer:.0f} bp (regularity {self.spacer_regularity:.2f}), "
            f"repeat identity {self.repeat_identity:.2f}, "
            f"spacer identity {self.spacer_identity:.2f}, "
            f"span {self.start}-{self.end}"
        )


def _seed_positions(dna: str, k: int, min_copies: int) -> dict[str, list[int]]:
    """k-mers occurring at least ``min_copies`` times."""
    pos: dict[str, list[int]] = defaultdict(list)
    for i in range(len(dna) - k + 1):
        pos[dna[i:i + k]].append(i)
    return {kmer: p for kmer, p in pos.items() if len(p) >= min_copies}


def _extend_unit(dna: str, positions: list[int], k: int, max_len: int,
                 max_mismatch_frac: float) -> tuple[str, int]:
    """Grow a seed k-mer rightwards while the copies stay similar."""
    length = k
    while length < max_len and positions[0] + length < len(dna):
        units = [dna[p:p + length + 1] for p in positions if p + length + 1 <= len(dna)]
        if len(units) < len(positions):
            break
        ref = units[0]
        if any(hamming(ref, u) / len(ref) > max_mismatch_frac for u in units[1:]):
            break
        length += 1
    return dna[positions[0]:positions[0] + length], length


def _consensus(units: list[str]) -> str:
    if not units:
        return ""
    width = min(len(u) for u in units)
    out: list[str] = []
    for j in range(width):
        counts: dict[str, int] = defaultdict(int)
        for u in units:
            counts[u[j]] += 1
        out.append(max(counts.items(), key=lambda kv: kv[1])[0])
    return "".join(out)


def _refine_positions(
    s: str,
    consensus: str,
    lo: int,
    hi: int,
    max_mismatch_frac: float,
) -> list[int]:
    """Re-scan a window for every approximate copy of ``consensus``.

    The seeding pass can only chain copies that share an exact k-mer, so a copy
    carrying a mutation inside the seed is skipped and the array is reported
    short -- or rejected outright, because the resulting double-length spacer
    falls outside the allowed range. This pass takes the provisional consensus
    back over the region and recovers the missed copies by Hamming distance,
    which is what CRISPR array finders do for the same reason.
    """
    L = len(consensus)
    if L == 0:
        return []
    lo = max(0, lo)
    hi = min(len(s), hi)
    if hi - lo < L:
        return []

    window = np.frombuffer(s[lo:hi].encode(), dtype=np.uint8)
    cons = np.frombuffer(consensus.encode(), dtype=np.uint8)
    n = len(window) - L + 1
    if n <= 0:
        return []
    # Sliding comparison: shape (n, L) view over the window, no copy.
    view = np.lib.stride_tricks.sliding_window_view(window, L)
    mismatches = (view != cons[None, :]).sum(axis=1)

    budget = int(max_mismatch_frac * L)
    order = np.argsort(mismatches, kind="stable")
    chosen: list[int] = []
    for off in order:
        if mismatches[off] > budget:
            break
        p = lo + int(off)
        if any(abs(p - q) < L for q in chosen):
            continue
        chosen.append(p)
    return sorted(chosen)


def _mean_pairwise_identity(seqs: list[str]) -> float:
    """Mean ungapped identity over all pairs, comparing the shared prefix."""
    if len(seqs) < 2:
        return 0.0
    scores: list[float] = []
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i], seqs[j]
            n = min(len(a), len(b))
            if n == 0:
                continue
            same = sum(1 for x in range(n) if a[x] == b[x])
            scores.append(same / max(len(a), len(b)))
    return mean(scores) if scores else 0.0


def find_repeat_arrays(
    dna: str,
    *,
    seed_k: int = 14,
    min_copies: int = 3,
    min_repeat: int = 18,
    max_repeat: int = 60,
    min_spacer: int = 15,
    max_spacer: int = 80,
    max_mismatch_frac: float = 0.15,
    both_strands: bool = True,
) -> list[RepeatArray]:
    """Find repeat-spacer arrays in a DNA sequence.

    Returns arrays sorted by copy number, highest first. Overlapping candidates
    are collapsed so that one physical array is reported once.
    """
    dna = dna.upper()
    results: list[RepeatArray] = []
    strands = [(1, dna)] + ([(-1, revcomp(dna))] if both_strands else [])

    for strand, s in strands:
        claimed: list[tuple[int, int]] = []
        seeds = _seed_positions(s, seed_k, min_copies)
        # Longest-recurring seeds first: they anchor the truest arrays.
        for kmer, positions in sorted(seeds.items(), key=lambda kv: -len(kv[1])):
            positions = sorted(positions)
            # Keep only positions whose gaps are inside the plausible period.
            chain = [positions[0]]
            for p in positions[1:]:
                gap = p - chain[-1]
                if min_repeat + min_spacer <= gap <= 2 * (max_repeat + max_spacer):
                    chain.append(p)
                elif gap > 2 * (max_repeat + max_spacer):
                    if len(chain) >= min_copies:
                        break
                    chain = [p]
            if len(chain) < min_copies:
                continue

            unit, ulen = _extend_unit(s, chain, seed_k, max_repeat, max_mismatch_frac)
            if not (min_repeat <= ulen <= max_repeat):
                continue

            cons = _consensus([s[p:p + ulen] for p in chain])

            # Recover copies the exact-seed pass missed, then keep the longest
            # run whose spacers all fall in the allowed range. Splitting on a
            # bad spacer rather than rejecting the whole candidate means a real
            # array flanked by a chance repeat is still reported.
            period = ulen + max_spacer
            refined = _refine_positions(
                s, cons,
                chain[0] - 2 * period, chain[-1] + 2 * period,
                max_mismatch_frac,
            ) or chain

            runs: list[list[int]] = [[refined[0]]]
            for prev, cur in zip(refined, refined[1:]):
                gap = cur - (prev + ulen)
                if min_spacer <= gap <= max_spacer:
                    runs[-1].append(cur)
                else:
                    runs.append([cur])
            run = max(runs, key=len)
            if len(run) < min_copies:
                continue

            units = [s[p:p + ulen] for p in run]
            cons = _consensus(units)
            spacers = [s[a + ulen:b] for a, b in zip(run, run[1:])]
            spacer_lens = [len(sp) for sp in spacers]
            chain = run

            start, end = chain[0], chain[-1] + ulen
            if any(start < ce and cs < end for cs, ce in claimed):
                continue
            claimed.append((start, end))

            rid = mean(
                1.0 - hamming(cons, u[:len(cons)]) / max(len(cons), 1) for u in units
            )
            arr = RepeatArray(
                start=start if strand == 1 else len(dna) - end,
                end=end if strand == 1 else len(dna) - start,
                repeat_consensus=cons,
                repeat_positions=chain,
                spacer_lengths=spacer_lens,
                spacers=spacers,
                repeat_identity=rid,
                spacer_identity=_mean_pairwise_identity(spacers),
                strand=strand,
            )
            results.append(arr)

    # One physical array can surface on both strands; keep the better call.
    results.sort(key=lambda a: (-a.copy_number, -a.repeat_identity))
    deduped: list[RepeatArray] = []
    for a in results:
        if any(a.start < b.end and b.start < a.end for b in deduped):
            continue
        deduped.append(a)
    return deduped


def array_transcript_units(array: RepeatArray) -> list[tuple[int, int]]:
    """Predicted short-RNA boundaries if the array is processed at repeat edges.

    The article notes the array is expressed as a set of distinct short RNAs.
    This returns the repeat-to-repeat intervals a processing enzyme cutting at
    repeat boundaries would release, which is the prediction a Northern blot or
    small-RNA sequencing run would test.
    """
    units: list[tuple[int, int]] = []
    rlen = array.repeat_length
    for a, b in zip(array.repeat_positions, array.repeat_positions[1:]):
        units.append((a, b + rlen))
    return units
