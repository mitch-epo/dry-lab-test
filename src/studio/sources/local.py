"""Source backed by a locally generated corpus.

This backend is what lets the studio run where the sequence databases are not
reachable. It answers exactly the same interface as the live backend, and it
refuses to hide what it is: its provenance block is marked synthetic and that
marking is propagated into every report the run produces.

Nothing in this module can see the planted ground truth. ``Corpus.read`` is
called without it, so a stage physically cannot consult the answers.
"""

from __future__ import annotations

import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import yaml

from ..bio.context import Gene, Locus, extract_locus
from ..bio.repeats import RepeatArray, find_repeat_arrays
from ..corpus.build import Corpus
from .base import ProteinRecord, Publication, SourceProvenance

_ASSETS = Path(__file__).resolve().parent.parent / "assets"


class LocalSource:
    """Serve sequences, loci and literature from a generated corpus."""

    def __init__(self, corpus_dir: Path | str, *, literature_path: Path | None = None):
        self.corpus_dir = Path(corpus_dir)
        # with_truth stays False: the planted answers are not loaded into the
        # object the stages hold.
        self.corpus = Corpus.read(self.corpus_dir, with_truth=False)
        self._literature_path = literature_path or (_ASSETS / "literature.yaml")
        self._array_cache: dict[str, list[RepeatArray]] = {}
        self._gene_index: dict[str, list] = {}
        for g in self.corpus.genes.values():
            self._gene_index.setdefault(g.contig_id, []).append(g)
        for contig in self._gene_index:
            self._gene_index[contig].sort(key=lambda g: g.start)

    # -- provenance ---------------------------------------------------------

    def provenance(self) -> SourceProvenance:
        return SourceProvenance(
            kind="local",
            description=f"generated corpus at {self.corpus_dir}",
            synthetic=bool(self.corpus.meta.get("synthetic", True)),
            corpus_meta={k: v for k, v in self.corpus.meta.items() if k != "references"},
            retrieved_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )

    # -- literature ---------------------------------------------------------

    def literature(self, query: str, *, limit: int = 20) -> list[Publication]:
        """Return the bundled digest entries relevant to ``query``.

        Matching is a keyword overlap against titles and claims. The digest
        carries no citations; see ``literature.yaml`` for why.
        """
        doc = yaml.safe_load(self._literature_path.read_text())
        terms = {t for t in query.lower().replace("_", " ").split() if len(t) > 2}
        out: list[tuple[int, Publication]] = []
        for entry in doc.get("entries", []):
            blob = " ".join(
                [entry.get("title", ""), *entry.get("claims", []), entry.get("observable", "")]
            ).lower()
            score = sum(1 for t in terms if t in blob)
            out.append((score, Publication(
                pub_id=entry["pub_id"],
                title=entry["title"],
                abstract=entry.get("observable", ""),
                claims=list(entry.get("claims", [])),
                source=f"offline-digest v{doc.get('digest_version')}",
            )))
        out.sort(key=lambda kv: -kv[0])
        return [p for _, p in out[:limit]]

    # -- sequences ----------------------------------------------------------

    def reference_proteins(self) -> dict[str, str]:
        return dict(self.corpus.meta.get("references", {}))

    def family_seed(self, family: str, *, limit: int = 12) -> list[str]:
        """Seed sequences for the family profile.

        For the local source the seed is drawn from the reference set plus a
        spread of corpus members of that family. Taking a spread rather than
        the most similar members is deliberate: a profile built from one tight
        clade will only ever find that clade.
        """
        refs = self.reference_proteins()
        seed: list[str] = []
        if family in refs:
            seed.append(refs[family])
        members = [g.protein for g in self.corpus.genes.values() if g.family == family]
        if members:
            members.sort(key=len)
            step = max(1, len(members) // max(1, limit - len(seed)))
            seed.extend(members[::step][: limit - len(seed)])
        return seed[:limit]

    def iter_proteins(self) -> Iterator[ProteinRecord]:
        for g in self.corpus.genes.values():
            contig = self.corpus.contigs.get(g.contig_id)
            yield ProteinRecord(
                protein_id=g.gene_id,
                sequence=g.protein,
                name=g.annotation or g.gene_id,
                lineage=contig.lineage if contig else "",
                contig_id=g.contig_id,
                source="local-corpus",
                annotation=g.annotation,
            )

    def protein(self, protein_id: str) -> ProteinRecord | None:
        g = self.corpus.genes.get(protein_id)
        if g is None:
            return None
        contig = self.corpus.contigs.get(g.contig_id)
        return ProteinRecord(
            protein_id=g.gene_id,
            sequence=g.protein,
            name=g.annotation or g.gene_id,
            lineage=contig.lineage if contig else "",
            contig_id=g.contig_id,
            source="local-corpus",
            annotation=g.annotation,
        )

    # -- genomic context ----------------------------------------------------

    def contig_dna(self, contig_id: str) -> str | None:
        c = self.corpus.contigs.get(contig_id)
        return c.dna if c else None

    def _arrays(self, contig_id: str) -> list[RepeatArray]:
        """Detect arrays on a contig, once, from the DNA.

        The detector runs over the sequence rather than reading the corpus's
        own array annotations: the pipeline must find the arrays the same way
        it would in real data, or the array evidence it reports would be
        circular.
        """
        if contig_id not in self._array_cache:
            dna = self.contig_dna(contig_id) or ""
            self._array_cache[contig_id] = find_repeat_arrays(dna) if dna else []
        return self._array_cache[contig_id]

    def _genes_as_bio(self, contig_id: str) -> list[Gene]:
        return [
            Gene(
                gene_id=g.gene_id, start=g.start, end=g.end, strand=g.strand,
                protein=g.protein, label="unk", annotation=g.annotation,
            )
            for g in self._gene_index.get(contig_id, [])
        ]

    def locus(self, protein_id: str, *, window_bp: int = 12_000) -> Locus | None:
        g = self.corpus.genes.get(protein_id)
        if g is None:
            return None
        contig = self.corpus.contigs.get(g.contig_id)
        return extract_locus(
            g.contig_id,
            self._genes_as_bio(g.contig_id),
            protein_id,
            window_bp=window_bp,
            arrays=self._arrays(g.contig_id),
            taxon=contig.lineage if contig else "",
            source="local-corpus",
        )

    def random_loci(self, n: int, *, seed: int = 0, window_bp: int = 12_000) -> list[Locus]:
        rng = random.Random(seed)
        ids = sorted(self.corpus.genes)
        picks = rng.sample(ids, min(n, len(ids)))
        out: list[Locus] = []
        for pid in picks:
            loc = self.locus(pid, window_bp=window_bp)
            if loc is not None:
                out.append(loc)
        return out
