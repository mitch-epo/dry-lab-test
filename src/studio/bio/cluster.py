"""Sequence clustering and phylogenetic spread.

Two numbers from this module decide whether a candidate is a real finding or a
sampling artefact:

* how many *non-redundant* sequences support it -- 400 near-identical proteins
  from one over-sequenced clade are one observation, not 400; and
* how far apart in the tree those sequences are -- a system recovered from
  distant hosts is under selection, while one confined to a single clade may be
  a recent, local accident.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from .align import align
from .search import reduce_alphabet


@dataclass
class Cluster:
    """A set of sequences represented by one member."""

    representative: str
    members: list[str] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.members)


def _kmer_set(seq: str, k: int = 4) -> set[str]:
    red = reduce_alphabet(seq)
    return {red[i:i + k] for i in range(len(red) - k + 1) if "X" not in red[i:i + k]}


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def greedy_cluster(
    records: dict[str, str],
    *,
    identity: float = 0.9,
    coverage: float = 0.8,
    prefilter_jaccard: float = 0.35,
) -> list[Cluster]:
    """Greedy incremental clustering, longest sequence first.

    Follows the CD-HIT/MMseqs2 shape: sort by length, and assign each sequence
    to the first representative it matches above the identity and coverage
    thresholds, otherwise start a new cluster. The k-mer Jaccard prefilter keeps
    the number of full alignments near-linear in practice.
    """
    order = sorted(records.items(), key=lambda kv: -len(kv[1]))
    clusters: list[Cluster] = []
    reps: list[tuple[str, str, set[str]]] = []  # (id, seq, kmers)

    for sid, seq in order:
        ks = _kmer_set(seq)
        placed = False
        for ci, (rid, rseq, rks) in enumerate(reps):
            if _jaccard(ks, rks) < prefilter_jaccard:
                continue
            a = align(seq, rseq)
            cov = a.aligned_len / max(1, min(len(seq), len(rseq)))
            if a.identity >= identity and cov >= coverage:
                clusters[ci].members.append(sid)
                placed = True
                break
        if not placed:
            clusters.append(Cluster(sid, [sid]))
            reps.append((sid, seq, ks))
    return clusters


def non_redundant(records: dict[str, str], *, identity: float = 0.9) -> list[str]:
    """Representative IDs after collapsing at ``identity``."""
    return [c.representative for c in greedy_cluster(records, identity=identity)]


@dataclass
class TaxonomicSpread:
    """How widely a candidate set is distributed across a taxonomy."""

    n_sequences: int
    n_clusters: int
    n_taxa: int
    n_phyla: int
    largest_clade_fraction: float
    lineages: dict[str, int] = field(default_factory=dict)

    @property
    def redundancy(self) -> float:
        """1.0 when every sequence is distinct; near 0 when they are one clade's copies."""
        return self.n_clusters / self.n_sequences if self.n_sequences else 0.0

    @property
    def is_clade_restricted(self) -> bool:
        return self.largest_clade_fraction > 0.8 or self.n_phyla <= 1

    def summary(self) -> str:
        return (
            f"{self.n_sequences} sequences -> {self.n_clusters} clusters at 90% id "
            f"across {self.n_taxa} taxa / {self.n_phyla} phyla; "
            f"largest clade holds {self.largest_clade_fraction:.0%}"
        )


def taxonomic_spread(
    records: dict[str, str],
    lineages: dict[str, str],
    *,
    identity: float = 0.9,
    phylum_depth: int = 1,
) -> TaxonomicSpread:
    """Combine sequence redundancy with taxonomic breadth.

    ``lineages`` maps sequence ID to a semicolon-delimited lineage string.
    """
    clusters = greedy_cluster(records, identity=identity)
    taxa = {lineages.get(sid, "unknown") for sid in records}
    phyla: dict[str, int] = defaultdict(int)
    for sid in records:
        parts = [p.strip() for p in lineages.get(sid, "unknown").split(";") if p.strip()]
        phyla[parts[phylum_depth] if len(parts) > phylum_depth else (parts[0] if parts else "unknown")] += 1
    largest = max(phyla.values()) / len(records) if records else 0.0
    return TaxonomicSpread(
        n_sequences=len(records),
        n_clusters=len(clusters),
        n_taxa=len(taxa),
        n_phyla=len(phyla),
        largest_clade_fraction=largest,
        lineages=dict(phyla),
    )
