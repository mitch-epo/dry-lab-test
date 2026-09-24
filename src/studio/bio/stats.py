"""Statistics used to defend, or demolish, a claim in a candidate report.

Every quantitative claim the workflow makes about a gene neighbourhood --
"this partner gene is specifically associated with the enzyme" -- is an
enrichment claim, and an enrichment claim without a null model is an anecdote.
These are the tests the reporting stage must cite and the triage stage checks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def _log_factorial(n: int) -> float:
    return math.lgamma(n + 1)


def _log_hypergeom(a: int, b: int, c: int, d: int) -> float:
    """log P of one 2x2 table under the hypergeometric null."""
    n = a + b + c + d
    return (
        _log_factorial(a + b) + _log_factorial(c + d)
        + _log_factorial(a + c) + _log_factorial(b + d)
        - _log_factorial(a) - _log_factorial(b)
        - _log_factorial(c) - _log_factorial(d) - _log_factorial(n)
    )


@dataclass
class FisherResult:
    """A 2x2 association test."""

    odds_ratio: float
    p_value: float
    table: tuple[int, int, int, int]
    n: int

    @property
    def significant(self) -> bool:
        return self.p_value < 0.05

    def summary(self) -> str:
        a, b, c, d = self.table
        return (
            f"OR={self.odds_ratio:.1f}, p={self.p_value:.2e} "
            f"(with/without: {a}/{b} in set, {c}/{d} in background)"
        )


def fisher_exact(a: int, b: int, c: int, d: int) -> FisherResult:
    """One-sided (enrichment) Fisher exact test on the table [[a, b], [c, d]].

    ``a`` = feature present in the candidate set, ``b`` = absent in the set,
    ``c`` = present in the background, ``d`` = absent in the background.
    Implemented directly so the studio does not depend on scipy.
    """
    n = a + b + c + d
    if n == 0:
        return FisherResult(1.0, 1.0, (a, b, c, d), 0)

    observed = _log_hypergeom(a, b, c, d)
    total = 0.0
    # Sum the tail: tables at least as extreme in the direction of enrichment.
    lo = max(0, a - d)
    hi = min(a + b, a + c)
    for k in range(a, hi + 1):
        kb = (a + b) - k
        kc = (a + c) - k
        kd = d - (k - a)
        if kb < 0 or kc < 0 or kd < 0:
            continue
        total += math.exp(_log_hypergeom(k, kb, kc, kd))
    del lo, observed

    num = (a + 0.5) * (d + 0.5)
    den = (b + 0.5) * (c + 0.5)
    return FisherResult(num / den, min(1.0, total), (a, b, c, d), n)


def benjamini_hochberg(pvals: list[float], alpha: float = 0.05) -> list[bool]:
    """Benjamini-Hochberg rejection flags at FDR ``alpha``.

    A campaign tests hundreds of neighbour labels for association. Reporting
    raw p-values across that many tests manufactures findings; the workflow
    corrects before any enrichment claim reaches a report.
    """
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    flags = [False] * n
    max_k = -1
    for rank, i in enumerate(order, start=1):
        if pvals[i] <= alpha * rank / n:
            max_k = rank
    if max_k > 0:
        for rank, i in enumerate(order, start=1):
            if rank <= max_k:
                flags[i] = True
    return flags


@dataclass
class PermutationResult:
    observed: float
    null_mean: float
    null_sd: float
    p_value: float
    n_permutations: int

    @property
    def z(self) -> float:
        return (self.observed - self.null_mean) / (self.null_sd or 1.0)

    def summary(self) -> str:
        return (
            f"observed {self.observed:.3f} vs null {self.null_mean:.3f}"
            f"+/-{self.null_sd:.3f} (z={self.z:.1f}, p={self.p_value:.3g}, "
            f"{self.n_permutations} permutations)"
        )


def permutation_test(
    observed: float,
    null_samples: list[float],
) -> PermutationResult:
    """Empirical p-value for ``observed`` against a sampled null distribution."""
    n = len(null_samples)
    if n == 0:
        return PermutationResult(observed, 0.0, 1.0, 1.0, 0)
    mu = sum(null_samples) / n
    var = sum((x - mu) ** 2 for x in null_samples) / n
    sd = math.sqrt(var)
    # Add-one correction so a p-value is never reported as exactly zero.
    extreme = sum(1 for x in null_samples if x >= observed)
    return PermutationResult(observed, mu, sd, (extreme + 1) / (n + 1), n)


def neighbour_enrichment(
    label: str,
    candidate_loci_labels: list[set[str]],
    background_loci_labels: list[set[str]],
) -> FisherResult:
    """Is ``label`` enriched beside the candidates relative to random loci?

    This is the test behind every "partner gene" claim. The background must be
    loci sampled from the same genomes, otherwise the test measures genome
    composition rather than association.
    """
    a = sum(1 for s in candidate_loci_labels if label in s)
    b = len(candidate_loci_labels) - a
    c = sum(1 for s in background_loci_labels if label in s)
    d = len(background_loci_labels) - c
    return fisher_exact(a, b, c, d)
