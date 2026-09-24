"""Pairwise alignment: BLOSUM62, Gotoh affine gaps, Karlin-Altschul statistics.

The dynamic program is vectorised along anti-diagonals with numpy. On an
anti-diagonal every cell depends only on the two preceding anti-diagonals, so
the whole Gotoh recurrence -- including the horizontal gap state, which is the
part that blocks naive row-wise vectorisation -- becomes three array
operations per diagonal.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .seq import AA as AA_LETTERS, AA_INDEX

# BLOSUM62 in the canonical NCBI residue order.
_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62_ROWS = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""


def _build_blosum62() -> np.ndarray:
    """21x21 score matrix indexed by :data:`studio.bio.seq.AA_INDEX`, with row/col
    20 reserved for ``X`` (scored as a mild mismatch)."""
    raw = np.array(
        [[int(v) for v in line.split()] for line in _B62_ROWS.strip().splitlines()],
        dtype=np.float32,
    )
    m = np.full((21, 21), -1.0, dtype=np.float32)
    for i, a in enumerate(_B62_ORDER):
        for j, b in enumerate(_B62_ORDER):
            m[AA_INDEX[a], AA_INDEX[b]] = raw[i, j]
    return m


BLOSUM62 = _build_blosum62()

# Karlin-Altschul parameters for BLOSUM62 with gap open 11 / extend 1, as
# published by NCBI BLAST. Used to convert raw scores to bits and E-values.
KA_LAMBDA = 0.267
KA_K = 0.041

GAP_OPEN = 11.0
GAP_EXTEND = 1.0

_NEG = np.float32(-1e9)


def encode(seq: str) -> np.ndarray:
    """Map a protein string to matrix indices; unknown residues become 20 (``X``)."""
    return np.fromiter((AA_INDEX.get(c, 20) for c in seq), dtype=np.int64, count=len(seq))


@dataclass
class Alignment:
    """Result of a pairwise alignment."""

    score: float
    bits: float
    evalue: float
    q_start: int
    q_end: int
    t_start: int
    t_end: int
    identity: float
    similarity: float
    gaps: int
    aligned_len: int
    q_aln: str = ""
    t_aln: str = ""

    @property
    def q_coverage_of(self) -> int:
        return self.q_end - self.q_start


def bit_score(raw: float) -> float:
    """Convert a raw alignment score to bits."""
    return (KA_LAMBDA * raw - np.log(KA_K)) / np.log(2.0)


def evalue(raw: float, m: int, n: int) -> float:
    """Karlin-Altschul E-value for a raw score in an m x n search space."""
    b = bit_score(raw)
    # Clamp the exponent: scores beyond ~1e3 bits underflow to 0 anyway.
    return float(m) * float(n) * float(2.0 ** -min(b, 1000.0))


def _fill(
    prof: np.ndarray,
    local: bool,
    gap_open: float,
    gap_extend: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fill the Gotoh H/E/F matrices along anti-diagonals.

    ``prof`` is the m x n matrix of position scores -- substitution scores for a
    sequence-sequence alignment, or PSSM column scores for a profile-sequence
    alignment. ``H`` is the best score of an alignment ending at (i, j) in any
    state, ``E`` the best ending with a gap in the query, ``F`` with a gap in
    the target.
    """
    m, n = prof.shape

    H = np.full((m + 1, n + 1), _NEG, dtype=np.float32)
    E = np.full((m + 1, n + 1), _NEG, dtype=np.float32)
    F = np.full((m + 1, n + 1), _NEG, dtype=np.float32)

    oe = np.float32(gap_open + gap_extend)
    ge = np.float32(gap_extend)

    if local:
        H[0, :] = 0.0
        H[:, 0] = 0.0
    else:
        # Semi-global: free gaps at the start of either sequence, so that a short
        # domain can align inside a long protein without paying terminal gaps.
        H[0, :] = 0.0
        H[:, 0] = 0.0

    prof = prof.astype(np.float32, copy=False)
    ii = np.arange(1, m + 1)
    for d in range(2, m + n + 1):
        jj = d - ii
        valid = (jj >= 1) & (jj <= n)
        if not valid.any():
            continue
        i_v = ii[valid]
        j_v = jj[valid]

        # Horizontal gap state: left neighbour (i, j-1) lives on diagonal d-1.
        e = np.maximum(H[i_v, j_v - 1] - oe, E[i_v, j_v - 1] - ge)
        # Vertical gap state: upper neighbour (i-1, j) also lives on diagonal d-1.
        f = np.maximum(H[i_v - 1, j_v] - oe, F[i_v - 1, j_v] - ge)
        # Match state: diagonal neighbour (i-1, j-1) lives on diagonal d-2.
        h = H[i_v - 1, j_v - 1] + prof[i_v - 1, j_v - 1]

        best = np.maximum(h, np.maximum(e, f))
        if local:
            best = np.maximum(best, 0.0)

        H[i_v, j_v] = best
        E[i_v, j_v] = e
        F[i_v, j_v] = f

    return H, E, F


def _traceback(
    qs: str,
    ts: str,
    H: np.ndarray,
    E: np.ndarray,
    F: np.ndarray,
    prof: np.ndarray,
    i: int,
    j: int,
    local: bool,
    gap_open: float,
    gap_extend: float,
) -> tuple[str, str, int, int]:
    """Walk back from (i, j) reconstructing the aligned strings."""
    oe = gap_open + gap_extend
    qa: list[str] = []
    ta: list[str] = []
    while i > 0 and j > 0:
        h = H[i, j]
        if local and h <= 0:
            break
        diag = H[i - 1, j - 1] + prof[i - 1, j - 1]
        if np.isclose(h, diag, atol=1e-3):
            qa.append(qs[i - 1])
            ta.append(ts[j - 1])
            i -= 1
            j -= 1
        elif np.isclose(h, E[i, j], atol=1e-3):
            # Unwind the horizontal gap run.
            while True:
                qa.append("-")
                ta.append(ts[j - 1])
                prev_open = H[i, j - 1] - oe
                prev_ext = E[i, j - 1] - gap_extend
                j -= 1
                if prev_open >= prev_ext or j == 0:
                    break
        else:
            while True:
                qa.append(qs[i - 1])
                ta.append("-")
                prev_open = H[i - 1, j] - oe
                prev_ext = F[i - 1, j] - gap_extend
                i -= 1
                if prev_open >= prev_ext or i == 0:
                    break
    return "".join(reversed(qa)), "".join(reversed(ta)), i, j


def align(
    query: str,
    target: str,
    *,
    local: bool = True,
    sub: np.ndarray | None = None,
    gap_open: float = GAP_OPEN,
    gap_extend: float = GAP_EXTEND,
    search_space: tuple[int, int] | None = None,
) -> Alignment:
    """Align two protein sequences with affine gaps.

    ``local`` selects Smith-Waterman; otherwise the alignment is semi-global
    (free terminal gaps), which is the right model for placing a domain inside a
    longer protein. ``search_space`` overrides the (m, n) used for the E-value so
    that a hit found in a database scan is scored against the whole database.
    """
    sub = BLOSUM62 if sub is None else sub
    q, t = encode(query), encode(target)
    if len(q) == 0 or len(t) == 0:
        return Alignment(0.0, 0.0, float("inf"), 0, 0, 0, 0, 0.0, 0.0, 0, 0)

    prof = sub[q[:, None], t[None, :]]
    H, E, F = _fill(prof, local, gap_open, gap_extend)

    if local:
        idx = int(np.argmax(H[1:, 1:]))
        i = idx // H[1:, 1:].shape[1] + 1
        j = idx % H[1:, 1:].shape[1] + 1
    else:
        # Best score on the last row or column: allows either sequence to hang off.
        last_row, last_col = H[len(q), 1:], H[1:, len(t)]
        if last_row.max() >= last_col.max():
            i, j = len(q), int(np.argmax(last_row)) + 1
        else:
            i, j = int(np.argmax(last_col)) + 1, len(t)

    raw = float(H[i, j])
    q_aln, t_aln, qi0, tj0 = _traceback(
        query, target, H, E, F, prof, i, j, local, gap_open, gap_extend
    )

    matches = sum(1 for a, b in zip(q_aln, t_aln) if a == b and a != "-")
    positives = sum(
        1 for a, b in zip(q_aln, t_aln)
        if a != "-" and b != "-" and sub[AA_INDEX.get(a, 20), AA_INDEX.get(b, 20)] > 0
    )
    gaps = sum(1 for a, b in zip(q_aln, t_aln) if a == "-" or b == "-")
    alen = len(q_aln) or 1

    m_sp, n_sp = search_space or (len(query), len(target))
    return Alignment(
        score=raw,
        bits=bit_score(raw),
        evalue=evalue(raw, m_sp, n_sp),
        q_start=qi0,
        q_end=i,
        t_start=tj0,
        t_end=j,
        identity=matches / alen,
        similarity=positives / alen,
        gaps=gaps,
        aligned_len=alen,
        q_aln=q_aln,
        t_aln=t_aln,
    )


def quick_score(query: str, target: str, *, sub: np.ndarray | None = None) -> float:
    """Local alignment score only, skipping traceback.

    Used where thousands of pairs are scored and only the ranking matters.
    """
    sub = BLOSUM62 if sub is None else sub
    q, t = encode(query), encode(target)
    if len(q) == 0 or len(t) == 0:
        return 0.0
    H, _, _ = _fill(sub[q[:, None], t[None, :]], True, GAP_OPEN, GAP_EXTEND)
    return float(H.max())


def profile_align(
    pssm: np.ndarray,
    target: str,
    *,
    local: bool = True,
    gap_open: float = GAP_OPEN,
    gap_extend: float = GAP_EXTEND,
    search_space: tuple[int, int] | None = None,
) -> Alignment:
    """Align a position-specific scoring matrix against a sequence.

    ``pssm`` has shape (L, 21) in :data:`studio.bio.seq.AA_INDEX` column order.
    The returned ``q_start``/``q_end`` are profile column coordinates, which is
    what the search stage uses to judge whether a hit spans the catalytic core
    of a domain or only clips its edge.
    """
    t = encode(target)
    if pssm.shape[0] == 0 or len(t) == 0:
        return Alignment(0.0, 0.0, float("inf"), 0, 0, 0, 0, 0.0, 0.0, 0, 0)

    prof = pssm[:, t]
    H, E, F = _fill(prof, local, gap_open, gap_extend)
    m = pssm.shape[0]

    if local:
        flat = H[1:, 1:]
        idx = int(np.argmax(flat))
        i, j = idx // flat.shape[1] + 1, idx % flat.shape[1] + 1
    else:
        last_row, last_col = H[m, 1:], H[1:, len(t)]
        if last_row.max() >= last_col.max():
            i, j = m, int(np.argmax(last_row)) + 1
        else:
            i, j = int(np.argmax(last_col)) + 1, len(t)

    raw = float(H[i, j])
    # Consensus of the profile stands in for the "query" string in traceback.
    consensus = "".join(AA_LETTERS[int(c)] for c in np.argmax(pssm[:, :20], axis=1))
    q_aln, t_aln, qi0, tj0 = _traceback(
        consensus, target, H, E, F, prof, i, j, local, gap_open, gap_extend
    )

    matches = sum(1 for a, b in zip(q_aln, t_aln) if a == b and a != "-")
    positives = sum(
        1 for k, (a, b) in enumerate(zip(q_aln, t_aln))
        if a != "-" and b != "-" and BLOSUM62[AA_INDEX.get(a, 20), AA_INDEX.get(b, 20)] > 0
    )
    gaps = sum(1 for a, b in zip(q_aln, t_aln) if a == "-" or b == "-")
    alen = len(q_aln) or 1
    m_sp, n_sp = search_space or (m, len(target))
    return Alignment(
        score=raw,
        bits=bit_score(raw),
        evalue=evalue(raw, m_sp, n_sp),
        q_start=qi0,
        q_end=i,
        t_start=tj0,
        t_end=j,
        identity=matches / alen,
        similarity=positives / alen,
        gaps=gaps,
        aligned_len=alen,
        q_aln=q_aln,
        t_aln=t_aln,
    )


def profile_scores_batch(
    pssm: np.ndarray,
    encoded: np.ndarray,
    lengths: np.ndarray,
    *,
    gap_open: float = GAP_OPEN,
    gap_extend: float = GAP_EXTEND,
) -> np.ndarray:
    """Best local alignment score of one profile against many sequences at once.

    The single-pair dynamic program spends almost all of its time in Python
    loop overhead rather than in arithmetic: an anti-diagonal of a 330x330
    problem is a handful of microseconds of numpy work wrapped in tens of
    microseconds of interpreter. Scoring a database one pair at a time pays
    that overhead once per pair.

    Batching inverts the cost. Every target advances through the same
    anti-diagonal at the same time, so the interpreter runs the loop once for
    the whole batch and each numpy operation covers ``T`` targets. The loop
    length is unchanged; the work per iteration grows; the overhead per target
    falls by roughly the batch size.

    Only the maximum score is returned -- no traceback -- which is all a
    database scan needs before deciding which pairs deserve a full alignment.
    Three rolling anti-diagonals are kept instead of full matrices, so memory
    is O(T * m) rather than O(T * m * n).

    ``encoded`` is (T, n_max) of residue indices padded with 20 (``X``);
    ``lengths`` is (T,) giving each target's true length.
    """
    pssm = pssm.astype(np.float32, copy=False)
    m = pssm.shape[0]
    T, n_max = encoded.shape
    if T == 0 or m == 0 or n_max == 0:
        return np.zeros(T, dtype=np.float32)

    oe = np.float32(gap_open + gap_extend)
    ge = np.float32(gap_extend)

    # Rolling buffers, indexed [target, profile position]; index 0 is the
    # boundary row and stays at zero for a local alignment.
    shape = (T, m + 1)
    h_prev1 = np.zeros(shape, dtype=np.float32)
    h_prev2 = np.zeros(shape, dtype=np.float32)
    e_prev = np.full(shape, _NEG, dtype=np.float32)
    f_prev = np.full(shape, _NEG, dtype=np.float32)
    best = np.zeros(T, dtype=np.float32)

    ii = np.arange(1, m + 1)
    lengths = lengths.reshape(T, 1)

    h_cur = np.zeros(shape, dtype=np.float32)
    e_cur = np.full(shape, _NEG, dtype=np.float32)
    f_cur = np.full(shape, _NEG, dtype=np.float32)

    for d in range(2, m + n_max + 1):
        jj = d - ii
        in_range = (jj >= 1) & (jj <= n_max)
        if not in_range.any():
            continue
        i_v = ii[in_range]
        j_v = jj[in_range]

        # A target shorter than j_v has no cell here.
        valid = j_v[None, :] <= lengths            # (T, k)

        # Profile column scores for this anti-diagonal, gathered per target.
        s = pssm[i_v[None, :] - 1, encoded[:, j_v - 1]]      # (T, k)

        e = np.maximum(h_prev1[:, i_v] - oe, e_prev[:, i_v] - ge)
        f = np.maximum(h_prev1[:, i_v - 1] - oe, f_prev[:, i_v - 1] - ge)
        h = h_prev2[:, i_v - 1] + s

        cell = np.maximum(np.maximum(h, e), f)
        np.maximum(cell, 0.0, out=cell)
        cell *= valid

        h_cur[:] = 0.0
        e_cur[:] = _NEG
        f_cur[:] = _NEG
        h_cur[:, i_v] = cell
        e_cur[:, i_v] = np.where(valid, e, _NEG)
        f_cur[:, i_v] = np.where(valid, f, _NEG)

        np.maximum(best, cell.max(axis=1), out=best)

        h_prev2, h_prev1, h_cur = h_prev1, h_cur, h_prev2
        e_prev, e_cur = e_cur, e_prev
        f_prev, f_cur = f_cur, f_prev

    return best


def encode_batch(seqs: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Pack sequences into a padded index matrix and a length vector."""
    if not seqs:
        return np.zeros((0, 0), dtype=np.int64), np.zeros(0, dtype=np.int64)
    n_max = max(len(s) for s in seqs)
    out = np.full((len(seqs), n_max), 20, dtype=np.int64)
    lengths = np.zeros(len(seqs), dtype=np.int64)
    for i, s in enumerate(seqs):
        e = encode(s)
        out[i, : len(e)] = e
        lengths[i] = len(e)
    return out, lengths
