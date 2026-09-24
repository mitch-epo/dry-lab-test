"""The data-source interface every stage reads through.

The workflow does not care whether a protein record came from UniProt over the
network or from a local snapshot, but it does care that the difference is
recorded. Every source reports a provenance block that is written into the run
manifest, so a campaign's conclusions can always be traced to the exact data
they rest on -- including the fact that a run used simulated data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Protocol, runtime_checkable

from ..bio.context import Gene, Locus


@dataclass
class Publication:
    """One piece of literature grounding."""

    pub_id: str
    title: str
    year: int | None = None
    abstract: str = ""
    claims: list[str] = field(default_factory=list)
    url: str = ""
    source: str = ""

    def cite(self) -> str:
        y = f", {self.year}" if self.year else ""
        return f"{self.title}{y} [{self.pub_id}]"


@dataclass
class ProteinRecord:
    """One protein plus the context needed to reason about it."""

    protein_id: str
    sequence: str
    name: str = ""
    lineage: str = ""
    contig_id: str = ""
    source: str = ""
    annotation: str = ""


@dataclass
class SourceProvenance:
    """Where a run's data came from. Written verbatim into the manifest."""

    kind: str                       # "local" or "live"
    description: str
    synthetic: bool
    endpoints: dict[str, str] = field(default_factory=dict)
    corpus_meta: dict = field(default_factory=dict)
    retrieved_at: str = ""

    def banner(self) -> str:
        if self.synthetic:
            return (
                "SIMULATED DATA -- every sequence in this run is generated. "
                "No biological claim in these outputs refers to a real organism."
            )
        return f"live data from {', '.join(self.endpoints) or self.description}"


@runtime_checkable
class DataSource(Protocol):
    """What the stages require of a source of sequences and literature."""

    def provenance(self) -> SourceProvenance: ...

    def literature(self, query: str, *, limit: int = 20) -> list[Publication]:
        """Literature grounding for a family or system."""
        ...

    def reference_proteins(self) -> dict[str, str]:
        """One representative protein per described family, for gene labelling."""
        ...

    def family_seed(self, family: str, *, limit: int = 12) -> list[str]:
        """Sequences to build the family profile from."""
        ...

    def iter_proteins(self) -> Iterator[ProteinRecord]:
        """Every protein in the searchable database."""
        ...

    def protein(self, protein_id: str) -> ProteinRecord | None: ...

    def locus(self, protein_id: str, *, window_bp: int = 12_000) -> Locus | None:
        """The genomic neighbourhood around a protein's gene, arrays included."""
        ...

    def contig_dna(self, contig_id: str) -> str | None: ...

    def random_loci(self, n: int, *, seed: int = 0) -> list[Locus]:
        """Background loci for enrichment tests.

        These must come from the same genomes as the candidates, otherwise an
        enrichment test measures genome composition instead of association.
        """
        ...


__all__ = [
    "DataSource", "Gene", "Locus", "ProteinRecord", "Publication", "SourceProvenance",
]
