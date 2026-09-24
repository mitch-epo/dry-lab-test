"""Source backed by the public sequence and literature databases.

Endpoints used, all public and unauthenticated:

* UniProtKB REST         -- protein records, sequences, taxonomy
* NCBI E-utilities       -- nucleotide records and genomic context
* InterPro / Pfam REST   -- family membership and domain boundaries
* Europe PMC REST        -- literature, with real identifiers

This backend is the reason the pipeline is written against an interface. The
studio's analysis stages are identical in both modes; only where the bytes come
from differs, and the run manifest records which was used.

NOTE ON AVAILABILITY. The environment this repository was developed in blocks
outbound connections to all four hosts at the egress proxy, so this backend is
written against the documented APIs but has not been executed end to end. The
request shapes, response field paths and pagination follow each service's
published REST contract. :func:`probe` reports reachability so a run fails
loudly at startup rather than silently halfway through, and
``studio doctor`` calls it.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Iterator

import httpx

from ..bio.context import Gene, Locus, extract_locus
from ..bio.repeats import find_repeat_arrays
from ..bio.seq import translate
from .base import ProteinRecord, Publication, SourceProvenance

UNIPROT = "https://rest.uniprot.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
INTERPRO = "https://www.ebi.ac.uk/interpro/api"
EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"

ENDPOINTS = {
    "uniprot": UNIPROT,
    "ncbi": EUTILS,
    "interpro": INTERPRO,
    "europepmc": EUROPEPMC,
}


@dataclass
class ProbeResult:
    endpoint: str
    url: str
    ok: bool
    status: int | None = None
    detail: str = ""


def probe(timeout: float = 10.0) -> list[ProbeResult]:
    """Check each endpoint. Run before a live campaign; ``studio doctor`` uses it."""
    checks = {
        "uniprot": f"{UNIPROT}/uniprotkb/P0A7B8.json?fields=accession",
        "ncbi": f"{EUTILS}/einfo.fcgi?db=protein&retmode=json",
        "interpro": f"{INTERPRO}/entry/pfam/PF00078/",
        "europepmc": f"{EUROPEPMC}/search?query=retron&format=json&pageSize=1",
    }
    out: list[ProbeResult] = []
    for name, url in checks.items():
        try:
            r = httpx.get(url, timeout=timeout, follow_redirects=True)
            out.append(ProbeResult(name, url, r.status_code == 200, r.status_code))
        except Exception as exc:  # network, TLS, proxy policy
            out.append(ProbeResult(name, url, False, None, f"{type(exc).__name__}: {exc}"))
    return out


class LiveSource:
    """Serve sequences, loci and literature from the public databases."""

    def __init__(
        self,
        *,
        family_query: str = "(xref:pfam-PF00078)",
        taxon: str | None = None,
        limit: int = 5000,
        email: str | None = None,
        timeout: float = 60.0,
        rate_delay: float = 0.34,   # NCBI asks for <=3 requests/second
        cache: dict | None = None,
    ):
        self.family_query = family_query
        self.taxon = taxon
        self.limit = limit
        self.email = email
        self.rate_delay = rate_delay
        self._client = httpx.Client(
            timeout=timeout,
            follow_redirects=True,
            headers={"User-Agent": "molecular-design-studio/0.1 (+dry-lab)"},
        )
        self._cache: dict[str, Any] = cache if cache is not None else {}
        self._proteins: dict[str, ProteinRecord] = {}
        self._last_call = 0.0

    # -- plumbing -----------------------------------------------------------

    def _get(self, url: str, params: dict | None = None) -> httpx.Response:
        delta = time.monotonic() - self._last_call
        if delta < self.rate_delay:
            time.sleep(self.rate_delay - delta)
        self._last_call = time.monotonic()
        r = self._client.get(url, params=params)
        r.raise_for_status()
        return r

    def close(self) -> None:
        self._client.close()

    def provenance(self) -> SourceProvenance:
        return SourceProvenance(
            kind="live",
            description=f"public databases; family query {self.family_query!r}",
            synthetic=False,
            endpoints=dict(ENDPOINTS),
            retrieved_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )

    # -- literature ---------------------------------------------------------

    def literature(self, query: str, *, limit: int = 20) -> list[Publication]:
        """Search Europe PMC. Returns records with real identifiers."""
        r = self._get(
            f"{EUROPEPMC}/search",
            {"query": query, "format": "json", "pageSize": min(limit, 100),
             "resultType": "core"},
        )
        out: list[Publication] = []
        for res in r.json().get("resultList", {}).get("result", []):
            pid = res.get("doi") or res.get("pmid") or res.get("id") or ""
            out.append(Publication(
                pub_id=f"doi:{pid}" if res.get("doi") else f"pmid:{pid}",
                title=res.get("title", ""),
                year=int(res["pubYear"]) if str(res.get("pubYear", "")).isdigit() else None,
                abstract=res.get("abstractText", "") or "",
                url=f"https://doi.org/{res['doi']}" if res.get("doi") else "",
                source="europepmc",
            ))
        return out

    # -- sequences ----------------------------------------------------------

    def _uniprot_stream(self, query: str, size: int, fields: str) -> Iterator[dict]:
        """Page through UniProt's stream endpoint following Link headers."""
        url = f"{UNIPROT}/uniprotkb/search"
        params = {"query": query, "format": "json", "size": min(size, 500), "fields": fields}
        fetched = 0
        while url and fetched < size:
            r = self._get(url, params)
            params = None  # subsequent pages carry their parameters in the Link URL
            payload = r.json()
            for rec in payload.get("results", []):
                yield rec
                fetched += 1
                if fetched >= size:
                    return
            link = r.headers.get("link", "")
            url = link.split(";")[0].strip("<> ") if 'rel="next"' in link else ""

    FIELDS = "accession,protein_name,sequence,organism_name,lineage,xref_embl,ft_domain"

    def _load_proteins(self) -> None:
        if self._proteins:
            return
        query = self.family_query
        if self.taxon:
            query = f"({query}) AND (taxonomy_id:{self.taxon})"
        for rec in self._uniprot_stream(query, self.limit, self.FIELDS):
            acc = rec.get("primaryAccession", "")
            seq = rec.get("sequence", {}).get("value", "")
            if not acc or not seq:
                continue
            org = rec.get("organism", {})
            lineage = ";".join(org.get("lineage", []) or []) or org.get("scientificName", "")
            name = (
                rec.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value", "")
            )
            embl = next(
                (x.get("id", "") for x in rec.get("uniProtKBCrossReferences", [])
                 if x.get("database") == "EMBL"),
                "",
            )
            self._proteins[acc] = ProteinRecord(
                protein_id=acc, sequence=seq, name=name, lineage=lineage,
                contig_id=embl, source="uniprot", annotation=name,
            )

    def reference_proteins(self) -> dict[str, str]:
        """One representative per described family, resolved through InterPro.

        Families are named by their InterPro or Pfam accession in live mode;
        the mapping below covers the RT-associated systems in the catalogue.
        """
        wanted = {
            "RT": "PF00078", "cas1": "PF01867", "cas2": "PF09827",
            "cas6": "PF01881", "avd": "PF06841", "retron_effector": "PF13884",
            "abi_partner": "PF07751", "maturase_x": "PF08388",
            "he_endonuclease": "PF01541", "integrase": "PF00589",
            "terminase": "PF03237", "portal": "PF04860", "transposase": "PF01609",
            "helicase": "PF00271", "methyltransferase": "PF01555",
        }
        out: dict[str, str] = {}
        for label, pfam in wanted.items():
            key = f"ref:{pfam}"
            if key in self._cache:
                out[label] = self._cache[key]
                continue
            try:
                r = self._get(
                    f"{UNIPROT}/uniprotkb/search",
                    {"query": f"(xref:pfam-{pfam}) AND (reviewed:true)",
                     "format": "json", "size": 1, "fields": "accession,sequence"},
                )
                results = r.json().get("results", [])
                if results:
                    seq = results[0].get("sequence", {}).get("value", "")
                    if seq:
                        self._cache[key] = seq
                        out[label] = seq
            except httpx.HTTPError:
                continue
        return out

    def family_seed(self, family: str, *, limit: int = 12) -> list[str]:
        """Reviewed members of the family, spread across the length distribution."""
        pfam = {"RT": "PF00078"}.get(family, family)
        q = f"(xref:pfam-{pfam}) AND (reviewed:true)" if pfam.startswith("PF") else family
        seqs = [
            rec.get("sequence", {}).get("value", "")
            for rec in self._uniprot_stream(q, limit * 4, "accession,sequence")
        ]
        seqs = [s for s in seqs if s]
        if not seqs:
            return []
        seqs.sort(key=len)
        step = max(1, len(seqs) // limit)
        return seqs[::step][:limit]

    def iter_proteins(self) -> Iterator[ProteinRecord]:
        self._load_proteins()
        yield from self._proteins.values()

    def protein(self, protein_id: str) -> ProteinRecord | None:
        self._load_proteins()
        if protein_id in self._proteins:
            return self._proteins[protein_id]
        try:
            r = self._get(f"{UNIPROT}/uniprotkb/{protein_id}.json")
        except httpx.HTTPError:
            return None
        rec = r.json()
        seq = rec.get("sequence", {}).get("value", "")
        if not seq:
            return None
        org = rec.get("organism", {})
        return ProteinRecord(
            protein_id=protein_id, sequence=seq,
            name=rec.get("proteinDescription", {}).get("recommendedName", {})
                    .get("fullName", {}).get("value", ""),
            lineage=";".join(org.get("lineage", []) or []),
            source="uniprot",
        )

    # -- genomic context ----------------------------------------------------

    def contig_dna(self, contig_id: str) -> str | None:
        """Fetch a nucleotide record from NCBI as FASTA."""
        key = f"dna:{contig_id}"
        if key in self._cache:
            return self._cache[key]
        params = {"db": "nuccore", "id": contig_id, "rettype": "fasta", "retmode": "text"}
        if self.email:
            params["email"] = self.email
        try:
            r = self._get(f"{EUTILS}/efetch.fcgi", params)
        except httpx.HTTPError:
            return None
        dna = "".join(line.strip() for line in r.text.splitlines() if not line.startswith(">"))
        self._cache[key] = dna
        return dna

    def _contig_features(self, contig_id: str) -> list[Gene]:
        """Pull CDS features for a nucleotide record from NCBI's JSON view."""
        key = f"feat:{contig_id}"
        if key in self._cache:
            return self._cache[key]
        params = {"db": "nuccore", "id": contig_id, "rettype": "ft", "retmode": "text"}
        if self.email:
            params["email"] = self.email
        try:
            r = self._get(f"{EUTILS}/efetch.fcgi", params)
        except httpx.HTTPError:
            return []

        # NCBI's feature-table format: a CDS line carries the coordinates, and
        # the qualifier lines that follow it carry product and protein_id.
        genes: list[Gene] = []
        for line in r.text.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2] == "CDS":
                a, b = parts[0].lstrip("<>"), parts[1].lstrip("<>")
                if not (a.isdigit() and b.isdigit()):
                    continue
                s, e = int(a), int(b)
                genes.append(Gene(
                    gene_id=f"{contig_id}:{min(s, e)}-{max(s, e)}",
                    start=min(s, e) - 1,
                    end=max(s, e),
                    strand=1 if s <= e else -1,
                ))
            elif genes and len(parts) >= 5 and parts[3] in ("product", "protein_id"):
                if parts[3] == "product":
                    genes[-1].annotation = parts[4]
                else:
                    genes[-1].gene_id = parts[4]

        dna = self.contig_dna(contig_id) or ""
        if dna:
            for g in genes:
                sub = dna[g.start:g.end]
                if g.strand == -1:
                    from ..bio.seq import revcomp
                    sub = revcomp(sub)
                g.protein = translate(sub)
        self._cache[key] = genes
        return genes

    def locus(self, protein_id: str, *, window_bp: int = 12_000) -> Locus | None:
        rec = self.protein(protein_id)
        if rec is None or not rec.contig_id:
            return None
        genes = self._contig_features(rec.contig_id)
        if not genes:
            return None
        anchor = min(
            genes,
            key=lambda g: 0 if g.gene_id == protein_id else
            (abs(len(g.protein) - len(rec.sequence)) + 10_000),
        )
        dna = self.contig_dna(rec.contig_id) or ""
        return extract_locus(
            rec.contig_id, genes, anchor.gene_id,
            window_bp=window_bp,
            arrays=find_repeat_arrays(dna) if dna else [],
            taxon=rec.lineage,
            source="ncbi",
        )

    def random_loci(self, n: int, *, seed: int = 0, window_bp: int = 12_000) -> list[Locus]:
        """Background loci drawn from the same contigs as the candidates."""
        import random as _random

        self._load_proteins()
        rng = _random.Random(seed)
        contigs = sorted({p.contig_id for p in self._proteins.values() if p.contig_id})
        out: list[Locus] = []
        for contig in rng.sample(contigs, min(len(contigs), n)):
            genes = self._contig_features(contig)
            if not genes:
                continue
            anchor = rng.choice(genes)
            dna = self.contig_dna(contig) or ""
            out.append(extract_locus(
                contig, genes, anchor.gene_id, window_bp=window_bp,
                arrays=find_repeat_arrays(dna) if dna else [], source="ncbi",
            ))
        return out
