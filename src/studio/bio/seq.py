"""Sequence primitives: alphabets, translation, composition, and low-complexity masking.

Nothing here talks to the network. Every downstream analysis stage is built on
these primitives so that the studio produces identical numbers whether it is
reading a live UniProt record or a local corpus snapshot.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}
NT = "ACGT"

# Standard genetic code (translation table 11 shares the amino-acid assignments;
# it differs only in which codons may act as alternative starts).
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn",
                            "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")

# Background amino-acid frequencies (Robinson & Robinson-style composition,
# renormalised over the 20-letter alphabet). Used as the null model for profile
# log-odds scores and for composition-bias flags.
BACKGROUND_AA = {
    "A": 0.0787, "C": 0.0157, "D": 0.0530, "E": 0.0660, "F": 0.0399,
    "G": 0.0693, "H": 0.0229, "I": 0.0593, "K": 0.0595, "L": 0.0960,
    "M": 0.0239, "N": 0.0414, "P": 0.0472, "Q": 0.0393, "R": 0.0534,
    "S": 0.0683, "T": 0.0541, "V": 0.0672, "W": 0.0120, "Y": 0.0305,
}


def revcomp(dna: str) -> str:
    """Reverse complement of a DNA string."""
    return dna.translate(_COMPLEMENT)[::-1]


def translate(dna: str, *, to_stop: bool = True) -> str:
    """Translate a DNA string in frame 0.

    Unknown or incomplete codons become ``X``; a trailing partial codon is
    dropped. With ``to_stop`` the protein is truncated at the first stop.
    """
    out: list[str] = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            if to_stop:
                break
            out.append("*")
            continue
        out.append(aa)
    return "".join(out)


def clean_protein(seq: str) -> str:
    """Uppercase a protein sequence and map non-standard letters to ``X``."""
    seq = seq.strip().upper()
    return "".join(c if c in AA_INDEX else "X" for c in seq)


def composition(seq: str) -> dict[str, float]:
    """Amino-acid frequencies over the standard alphabet (``X`` excluded)."""
    counts = Counter(c for c in seq if c in AA_INDEX)
    total = sum(counts.values()) or 1
    return {a: counts.get(a, 0) / total for a in AA}


def shannon_entropy(seq: str) -> float:
    """Shannon entropy in bits over the residues present in ``seq``."""
    counts = Counter(c for c in seq if c in AA_INDEX)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counts.values())


def low_complexity_fraction(seq: str, window: int = 20, threshold: float = 3.0) -> float:
    """Fraction of ``window``-residue windows whose entropy falls below ``threshold``.

    A SEG-style heuristic. High values mark sequences whose profile hits are
    likely to be composition artefacts rather than genuine homology, and the
    triage stage treats that as a kill criterion.
    """
    if len(seq) < window:
        return 1.0 if shannon_entropy(seq) < threshold else 0.0
    flagged = sum(
        1 for i in range(len(seq) - window + 1)
        if shannon_entropy(seq[i:i + window]) < threshold
    )
    return flagged / (len(seq) - window + 1)


def gc_content(dna: str) -> float:
    """GC fraction of a DNA string."""
    dna = dna.upper()
    n = sum(1 for c in dna if c in "ACGT")
    if n == 0:
        return 0.0
    return sum(1 for c in dna if c in "GC") / n


def hamming(a: str, b: str) -> int:
    """Hamming distance; sequences of unequal length are compared over the overlap
    and the length difference is added."""
    d = sum(1 for x, y in zip(a, b) if x != y)
    return d + abs(len(a) - len(b))


def identity(a: str, b: str) -> float:
    """Ungapped identity over the overlap of two equal-ish length strings."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / max(len(a), len(b))


@dataclass(frozen=True)
class Feature:
    """A located sub-sequence: a CDS, a repeat unit, a spacer, a domain call."""

    kind: str
    start: int          # 0-based, inclusive
    end: int            # 0-based, exclusive
    strand: int = 1     # +1 or -1
    label: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlaps(self, other: "Feature") -> bool:
        return self.start < other.end and other.start < self.end


def parse_fasta(text: str) -> dict[str, str]:
    """Minimal FASTA parser returning ``{id: sequence}`` keyed on the first token."""
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name = line[1:].strip().split()[0] if line[1:].strip() else ""
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def write_fasta(records: dict[str, str], width: int = 60) -> str:
    """Serialise ``{id: sequence}`` to FASTA text."""
    lines: list[str] = []
    for name, seq in records.items():
        lines.append(f">{name}")
        for i in range(0, len(seq), width):
            lines.append(seq[i:i + width])
    return "\n".join(lines) + "\n"


_ORF_RE = re.compile(r"(?=(ATG(?:...)*?(?:TAA|TAG|TGA)))")


def find_orfs(dna: str, min_aa: int = 50) -> list[Feature]:
    """Find ORFs on both strands, returning them in forward-strand coordinates.

    Used when the studio is handed raw contigs with no annotation.
    """
    found: list[Feature] = []
    n = len(dna)
    for strand in (1, -1):
        s = dna if strand == 1 else revcomp(dna)
        for m in _ORF_RE.finditer(s.upper()):
            orf = m.group(1)
            if len(orf) // 3 - 1 < min_aa:
                continue
            a, b = m.start(1), m.start(1) + len(orf)
            if strand == 1:
                found.append(Feature("CDS", a, b, 1))
            else:
                found.append(Feature("CDS", n - b, n - a, -1))
    return sorted(found, key=lambda f: (f.start, f.end))
