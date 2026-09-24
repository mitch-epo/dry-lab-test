"""Database scanning: k-mer prefilter, ungapped extension, then full alignment.

Full Gotoh alignment against a whole database is not affordable -- the campaign
in the article screened on the order of 10^5 proteins. This module uses the
same three-stage funnel BLAST and MMseqs2 use:

1. a reduced-alphabet k-mer index proposes targets and the alignment diagonals
   the shared words fall on;
2. an ungapped extension along those diagonals -- one array gather each --
   discards the large majority of proposals; and
3. only what survives pays for the quadratic gapped alignment.

Requiring shared k-mers to agree on a diagonal is what suppresses the
composition-driven matches a raw word count ranks highly.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from .align import BLOSUM62, Alignment, align, encode, profile_scores_batch
from .profile import Profile
from .seq import clean_protein

# Murphy et al. 10-letter reduction: groups residues that interchange freely,
# which is what lets the prefilter see homology below ~30% identity.
MURPHY10 = {
    **{c: "L" for c in "LVIM"},
    "C": "C",
    "A": "A",
    "G": "G",
    **{c: "S" for c in "ST"},
    "P": "P",
    **{c: "F" for c in "FYW"},
    **{c: "E" for c in "EDNQ"},
    **{c: "K" for c in "KR"},
    "H": "H",
}


def reduce_alphabet(seq: str) -> str:
    """Project a protein onto the Murphy-10 alphabet; unknown residues become ``X``."""
    return "".join(MURPHY10.get(c, "X") for c in seq)


@dataclass
class Hit:
    """One scored database hit."""

    target_id: str
    alignment: Alignment
    prefilter_score: int
    ungapped_score: float = 0.0
    z: float = 0.0
    core_coverage: float = 1.0

    @property
    def bits(self) -> float:
        return self.alignment.bits

    @property
    def evalue(self) -> float:
        return self.alignment.evalue


@dataclass
class ScanStats:
    """What the funnel discarded at each stage -- reported so a campaign can
    show that a small hit list reflects real filtering, not a broken index."""

    database: int = 0
    prefilter_passed: int = 0
    ungapped_passed: int = 0
    aligned: int = 0
    reported: int = 0

    def summary(self) -> str:
        return (
            f"{self.database} in database -> {self.prefilter_passed} past k-mer "
            f"prefilter -> {self.ungapped_passed} past ungapped extension -> "
            f"{self.aligned} aligned -> {self.reported} reported"
        )


class SequenceIndex:
    """Inverted k-mer index over a protein database."""

    def __init__(self, k: int = 4) -> None:
        self.k = k
        self.ids: list[str] = []
        self.seqs: list[str] = []
        self._enc: list[np.ndarray] = []
        self._pos: dict[str, list[tuple[int, int]]] = defaultdict(list)
        self._total_residues = 0
        self.last_scan = ScanStats()

    def add(self, seq_id: str, seq: str) -> None:
        seq = clean_protein(seq)
        idx = len(self.ids)
        self.ids.append(seq_id)
        self.seqs.append(seq)
        self._enc.append(encode(seq))
        self._total_residues += len(seq)
        red = reduce_alphabet(seq)
        for i in range(len(red) - self.k + 1):
            kmer = red[i:i + self.k]
            if "X" in kmer:
                continue
            self._pos[kmer].append((idx, i))

    def add_many(self, records: dict[str, str]) -> None:
        for sid, s in records.items():
            self.add(sid, s)

    def __len__(self) -> int:
        return len(self.ids)

    def get(self, seq_id: str) -> str:
        return self.seqs[self.ids.index(seq_id)]

    @property
    def total_residues(self) -> int:
        return self._total_residues

    def search_space(self, query_length: int) -> tuple[int, int]:
        """Karlin-Altschul search space for a query of this length.

        E = m*n*2^-S' takes m as the QUERY length and n as the database size in
        residues. Substituting the database's sequence *count* for m -- an easy
        slip, since both are properties of the database -- makes the E-value
        independent of query length, which removes the very length correction
        the statistic exists to provide and systematically mis-scores short
        proteins relative to long ones.
        """
        return (max(1, query_length), self._total_residues)

    # -- stage 1: k-mer prefilter -------------------------------------------

    def _prefilter(
        self, query: str, top: int, n_diagonals: int = 3
    ) -> list[tuple[int, int, list[int]]]:
        """Return ``(target_index, best k-mer count, top diagonals)`` per target.

        Keeping the diagonal offsets, not just the count, is what lets the
        caller run a cheap ungapped extension before committing to the full
        dynamic program.
        """
        red = reduce_alphabet(clean_protein(query))
        diag: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        for qi in range(len(red) - self.k + 1):
            kmer = red[qi:qi + self.k]
            if "X" in kmer:
                continue
            for tidx, ti in self._pos.get(kmer, ()):
                diag[tidx][ti - qi] += 1

        scored: list[tuple[int, int, list[int]]] = []
        for tidx, d in diag.items():
            best = sorted(d.items(), key=lambda kv: -kv[1])[:n_diagonals]
            scored.append((tidx, best[0][1], [off for off, _ in best]))
        scored.sort(key=lambda x: -x[1])
        return scored[:top]

    # -- stage 2: ungapped extension ----------------------------------------

    @staticmethod
    def _max_segment(scores: np.ndarray) -> float:
        """Maximum-scoring contiguous segment, vectorised.

        ``max_j (cumsum[j] - min_{i<=j} cumsum[i])`` is Kadane's algorithm
        written as two accumulate calls.
        """
        if scores.size == 0:
            return 0.0
        c = np.cumsum(scores, dtype=np.float64)
        prefix_min = np.minimum.accumulate(np.concatenate(([0.0], c[:-1])))
        return float(np.max(c - prefix_min))

    def _ungapped_best(
        self, matrix: np.ndarray, t_idx: int, diagonals: list[int], q_len: int
    ) -> float:
        """Best ungapped segment score for a target over the given diagonals.

        ``matrix`` has shape (q_len, 21): PSSM rows for a profile search, or
        substitution-matrix rows gathered for the query in a sequence search.
        """
        t = self._enc[t_idx]
        best = 0.0
        for off in diagonals:
            i_lo = max(0, -off)
            i_hi = min(q_len, len(t) - off)
            if i_hi <= i_lo:
                continue
            qi = np.arange(i_lo, i_hi)
            best = max(best, self._max_segment(matrix[qi, t[qi + off]]))
        return best

    def _batch_scores(self, pssm: np.ndarray, indices: list[int], batch: int = 256) -> np.ndarray:
        """Best local score of a profile against the given targets, batched.

        Targets are grouped by length so that padding a batch to its longest
        member costs little; an unsorted batch would pad every short sequence
        out to the longest in the database and throw the saving away.
        """
        if not indices:
            return np.zeros(0, dtype=np.float32)
        order = sorted(range(len(indices)), key=lambda k: len(self.seqs[indices[k]]))
        out = np.zeros(len(indices), dtype=np.float32)
        for start in range(0, len(order), batch):
            chunk = order[start:start + batch]
            enc, lengths = self._padded([indices[k] for k in chunk])
            scores = profile_scores_batch(pssm, enc, lengths)
            for k, sc in zip(chunk, scores):
                out[k] = sc
        return out

    def _padded(self, indices: list[int]) -> tuple[np.ndarray, np.ndarray]:
        n_max = max(len(self._enc[i]) for i in indices)
        enc = np.full((len(indices), n_max), 20, dtype=np.int64)
        lengths = np.zeros(len(indices), dtype=np.int64)
        for row, i in enumerate(indices):
            e = self._enc[i]
            enc[row, : len(e)] = e
            lengths[row] = len(e)
        return enc, lengths

    # -- stage 3: full alignment --------------------------------------------

    def search(
        self,
        query: str,
        *,
        max_evalue: float = 1e-3,
        prefilter_top: int = 2000,
        min_prefilter: int = 2,
        ungapped_cutoff: float | None = None,
        limit: int = 200,
    ) -> list[Hit]:
        """Sequence-vs-database search.

        ``ungapped_cutoff`` defaults to the raw score corresponding to roughly
        11 bits, below which no alignment can reach a useful E-value in a
        database of this size.
        """
        query = clean_protein(query)
        q = encode(query)
        rows = BLOSUM62[q]  # (len(query), 21)
        space = self.search_space(len(query))
        if ungapped_cutoff is None:
            ungapped_cutoff = 30.0

        stats = ScanStats(database=len(self.ids))
        survivors: list[tuple[int, int, float]] = []
        for tidx, pre, diagonals in self._prefilter(query, prefilter_top):
            if pre < min_prefilter:
                continue
            stats.prefilter_passed += 1
            ung = self._ungapped_best(rows, tidx, diagonals, len(query))
            if ung < ungapped_cutoff:
                continue
            stats.ungapped_passed += 1
            survivors.append((tidx, pre, ung))

        # Score every survivor in one batched pass, then pay for a traceback
        # only where the score could clear the threshold.
        hits: list[Hit] = []
        if survivors:
            scores = self._batch_scores(rows, [t for t, _, _ in survivors])
            from .align import bit_score, evalue as _evalue
            for (tidx, pre, ung), raw in zip(survivors, scores):
                if _evalue(float(raw), *space) > max_evalue:
                    continue
                a = align(query, self.seqs[tidx], search_space=space)
                stats.aligned += 1
                if a.evalue <= max_evalue:
                    hits.append(Hit(self.ids[tidx], a, pre, ungapped_score=ung))
        hits.sort(key=lambda h: h.alignment.score, reverse=True)
        stats.reported = min(len(hits), limit)
        self.last_scan = stats
        return hits[:limit]

    def profile_search(
        self,
        profile: Profile,
        *,
        max_evalue: float = 1e-3,
        min_z: float = 6.0,
        min_core_coverage: float = 0.7,
        prefilter_top: int = 20_000,
        min_prefilter: int = 2,
        ungapped_fraction: float = 0.5,
        limit: int = 5000,
    ) -> list[Hit]:
        """Profile-vs-database search with an explicit core-coverage requirement.

        A hit is kept only if it clears the E-value threshold, sits ``min_z``
        standard deviations above the profile's shuffled-decoy null, *and*
        covers at least ``min_core_coverage`` of the profile's high-information
        columns. The last condition separates a real domain call from a
        fragment matching only the domain's loosely conserved flanks.

        The ungapped gate is set to ``ungapped_fraction`` of the raw score that
        ``min_z`` demands. Gapped alignment can only improve on the best
        ungapped diagonal, so a target far below that bound cannot reach the
        threshold and is dropped without a full alignment.
        """
        space = self.search_space(profile.length)
        core = set(profile.core_columns)
        gate = ungapped_fraction * (profile.decoy_mu + min_z * profile.decoy_sigma)

        stats = ScanStats(database=len(self.ids))
        survivors: list[tuple[int, int, float]] = []
        for tidx, pre, diagonals in self._prefilter(profile.consensus, prefilter_top):
            if pre < min_prefilter:
                continue
            stats.prefilter_passed += 1
            ung = self._ungapped_best(profile.pssm, tidx, diagonals, profile.length)
            if ung < gate:
                continue
            stats.ungapped_passed += 1
            survivors.append((tidx, pre, ung))

        # One batched scoring pass over every survivor. The full alignment --
        # needed for the coordinates that core coverage is computed from -- is
        # run only where the batched score already clears both thresholds.
        hits: list[Hit] = []
        if survivors:
            from .align import evalue as _evalue

            scores = self._batch_scores(profile.pssm, [t for t, _, _ in survivors])
            for (tidx, pre, ung), raw in zip(survivors, scores):
                raw = float(raw)
                if profile.z(raw) < min_z or _evalue(raw, *space) > max_evalue:
                    continue
                a = profile.search(self.seqs[tidx], search_space=space)
                stats.aligned += 1
                if a.evalue > max_evalue:
                    continue
                z = profile.z(a.score)
                if z < min_z:
                    continue
                cov = (
                    len([c for c in core if a.q_start <= c < a.q_end]) / len(core)
                    if core else 1.0
                )
                if cov < min_core_coverage:
                    continue
                hits.append(
                    Hit(self.ids[tidx], a, pre, ungapped_score=ung, z=z, core_coverage=cov)
                )
        hits.sort(key=lambda h: h.alignment.score, reverse=True)
        stats.reported = min(len(hits), limit)
        self.last_scan = stats
        return hits[:limit]


    def profile_scan(
        self,
        profile: Profile,
        *,
        min_z: float = 5.0,
        prefilter_top: int = 20_000,
        min_prefilter: int = 2,
        ungapped_fraction: float = 0.5,
    ) -> list[tuple[str, float, float]]:
        """Score a profile against the database without any traceback.

        Returns ``(target_id, z, raw_score)`` for every target clearing
        ``min_z``. This is the funnel's first three stages only: k-mer
        prefilter, ungapped extension, batched scoring. Callers that need
        alignment coordinates -- and so a traceback -- run
        :func:`studio.bio.align.profile_align` on the few targets they
        actually care about, rather than paying for one per pair scanned.
        """
        gate = ungapped_fraction * (profile.decoy_mu + min_z * profile.decoy_sigma)
        survivors: list[int] = []
        for tidx, pre, diagonals in self._prefilter(profile.consensus, prefilter_top):
            if pre < min_prefilter:
                continue
            if self._ungapped_best(profile.pssm, tidx, diagonals, profile.length) < gate:
                continue
            survivors.append(tidx)
        if not survivors:
            return []
        scores = self._batch_scores(profile.pssm, survivors)
        out: list[tuple[str, float, float]] = []
        for tidx, raw in zip(survivors, scores):
            z = profile.z(float(raw))
            if z >= min_z:
                out.append((self.ids[tidx], z, float(raw)))
        out.sort(key=lambda t: -t[1])
        return out


def reciprocal_best_hits(
    index_a: SequenceIndex,
    index_b: SequenceIndex,
    *,
    max_evalue: float = 1e-5,
) -> dict[str, str]:
    """Reciprocal best hits between two databases.

    The workflow uses this to ask whether a candidate's partner gene is simply a
    diverged copy of a protein already present in a described system -- the most
    common way an apparently novel architecture turns out to be a known one.
    """
    best_ab: dict[str, tuple[str, float]] = {}
    for sid, seq in zip(index_a.ids, index_a.seqs):
        hs = index_b.search(seq, max_evalue=max_evalue, limit=1)
        if hs:
            best_ab[sid] = (hs[0].target_id, hs[0].alignment.score)

    rbh: dict[str, str] = {}
    for sid, (tid, _) in best_ab.items():
        back = index_a.search(index_b.get(tid), max_evalue=max_evalue, limit=1)
        if back and back[0].target_id == sid:
            rbh[sid] = tid
    return rbh
