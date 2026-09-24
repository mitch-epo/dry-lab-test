"""Genomic neighbourhood analysis and system-architecture comparison.

"Family members or genomic neighbours that fit no described system" is the
search criterion in the article, and it is a statement about architecture, not
about sequence. This module turns a locus into a canonical architecture string
-- ordered, strand-aware gene labels plus non-coding components -- and scores
that string against a catalogue of described systems.

A candidate becomes interesting when it is confidently inside the family and
its architecture matches nothing in the catalogue. Both halves matter: without
the first the hit is noise, and without the second it is a rediscovery.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import mean

from .repeats import RepeatArray


@dataclass
class Gene:
    """One annotated CDS on a contig."""

    gene_id: str
    start: int
    end: int
    strand: int
    protein: str = ""
    label: str = "unk"           # functional label, "unk" when uncharacterised
    annotation: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def aa_length(self) -> int:
        return len(self.protein)


@dataclass
class Locus:
    """A contig region centred on a gene of interest."""

    contig_id: str
    anchor: Gene
    genes: list[Gene] = field(default_factory=list)
    arrays: list[RepeatArray] = field(default_factory=list)
    taxon: str = ""
    source: str = ""

    @property
    def start(self) -> int:
        return min([g.start for g in self.genes] + [a.start for a in self.arrays])

    @property
    def end(self) -> int:
        return max([g.end for g in self.genes] + [a.end for a in self.arrays])


@dataclass
class Neighbour:
    gene: Gene
    offset_genes: int         # signed position relative to the anchor
    intergenic_gap: int       # bp between this gene and the one before it
    same_strand: bool
    co_directional_run: bool  # part of an uninterrupted same-strand run with the anchor


def extract_locus(
    contig_id: str,
    genes: list[Gene],
    anchor_id: str,
    *,
    window_bp: int = 12_000,
    max_genes_each_side: int = 6,
    arrays: list[RepeatArray] | None = None,
    taxon: str = "",
    source: str = "",
) -> Locus:
    """Pull the neighbourhood around ``anchor_id`` out of a contig's gene list."""
    ordered = sorted(genes, key=lambda g: g.start)
    idx = next((i for i, g in enumerate(ordered) if g.gene_id == anchor_id), None)
    if idx is None:
        raise KeyError(f"{anchor_id} not on contig {contig_id}")
    anchor = ordered[idx]

    lo = idx
    while lo > 0 and idx - lo < max_genes_each_side \
            and anchor.start - ordered[lo - 1].end <= window_bp:
        lo -= 1
    hi = idx
    while hi < len(ordered) - 1 and hi - idx < max_genes_each_side \
            and ordered[hi + 1].start - anchor.end <= window_bp:
        hi += 1

    window = ordered[lo:hi + 1]
    span_start = min(g.start for g in window) - 500
    span_end = max(g.end for g in window) + 500
    kept_arrays = [
        a for a in (arrays or []) if a.start < span_end and span_start < a.end
    ]
    return Locus(contig_id, anchor, window, kept_arrays, taxon=taxon, source=source)


def neighbours(locus: Locus) -> list[Neighbour]:
    """Describe each gene in the locus relative to the anchor."""
    ordered = sorted(locus.genes, key=lambda g: g.start)
    a_idx = next(i for i, g in enumerate(ordered) if g.gene_id == locus.anchor.gene_id)
    out: list[Neighbour] = []
    for i, g in enumerate(ordered):
        gap = g.start - ordered[i - 1].end if i > 0 else 0
        same = g.strand == locus.anchor.strand
        # Co-directional run: every gene between here and the anchor shares the
        # anchor's strand, the usual operon proxy when no transcriptomics exist.
        span = ordered[min(i, a_idx):max(i, a_idx) + 1]
        run = all(x.strand == locus.anchor.strand for x in span)
        out.append(Neighbour(g, i - a_idx, gap, same, run))
    return out


def predicted_operon(locus: Locus, max_gap: int = 60) -> list[Gene]:
    """Genes plausibly co-transcribed with the anchor.

    Same strand, uninterrupted, and separated by no more than ``max_gap`` bp.
    60 bp is the conventional bacterial cut-off; genes further apart are usually
    separately transcribed.
    """
    ordered = sorted(locus.genes, key=lambda g: g.start)
    a_idx = next(i for i, g in enumerate(ordered) if g.gene_id == locus.anchor.gene_id)
    keep = [ordered[a_idx]]

    i = a_idx
    while i > 0:
        prev, cur = ordered[i - 1], ordered[i]
        if prev.strand != locus.anchor.strand or cur.start - prev.end > max_gap:
            break
        keep.insert(0, prev)
        i -= 1
    i = a_idx
    while i < len(ordered) - 1:
        cur, nxt = ordered[i], ordered[i + 1]
        if nxt.strand != locus.anchor.strand or nxt.start - cur.end > max_gap:
            break
        keep.append(nxt)
        i += 1
    return keep


def architecture(
    locus: Locus,
    *,
    include_arrays: bool = True,
    array_detail: str = "full",
) -> str:
    """Canonical architecture string for a locus.

    Genes appear in coordinate order as ``label(+/-)``, with repeat arrays
    inserted at their position. Strands are expressed relative to the anchor
    so that a locus and its reverse complement produce the same string.

    ``array_detail`` selects how much of an array is written into the string:
    ``"full"`` gives its geometry (for display), ``"class"`` gives only its
    architectural class (for grouping). Grouping on geometry is a mistake with
    teeth -- copy number and repeat length differ between instances of the
    same system, so every locus becomes its own group of one and then fails
    every check that asks whether the architecture recurs. The component that
    makes the system most recognisable ends up being what disqualifies it.
    """
    flip = locus.anchor.strand == -1
    items: list[tuple[int, str]] = []
    for g in locus.genes:
        rel = g.strand * (-1 if flip else 1)
        items.append((g.start, f"{g.label}({'+' if rel > 0 else '-'})"))
    if include_arrays:
        for a in locus.arrays:
            token = (
                f"array[{a.copy_number}x{a.repeat_length}]" if array_detail == "full"
                else f"array[{a.class_token()}]"
            )
            items.append((a.start, token))

    items.sort(key=lambda kv: kv[0], reverse=flip)
    tokens = [text for _, text in items]

    if array_detail == "class":
        # The detector sometimes splits one physical array into two calls, or
        # reports a short spurious neighbour beside a real one. At class level
        # that is noise, not architecture, so adjacent identical array tokens
        # collapse -- otherwise one locus lands in a group of its own and then
        # fails every check that asks whether its architecture recurs.
        collapsed: list[str] = []
        for t in tokens:
            if t.startswith("array[") and collapsed and collapsed[-1] == t:
                continue
            collapsed.append(t)
        tokens = collapsed

    return "-".join(tokens)


def core_architecture(
    locus: Locus, *, max_gap: int = 60, array_detail: str = "full"
) -> str:
    """Architecture restricted to the predicted operon plus adjacent arrays.

    Comparing full windows to catalogue entries produces spurious mismatches,
    because the window includes housekeeping genes that have nothing to do with
    the system. The operon-scoped string is the fair unit of comparison.
    """
    op = predicted_operon(locus, max_gap=max_gap)
    sub = Locus(locus.contig_id, locus.anchor, op, locus.arrays, locus.taxon, locus.source)
    return architecture(sub, array_detail=array_detail)


@dataclass
class KnownSystem:
    """A described system, as the catalogue records it."""

    system_id: str
    name: str
    required_labels: list[str]
    forbidden_labels: list[str] = field(default_factory=list)
    optional_labels: list[str] = field(default_factory=list)
    requires_array: bool = False
    forbids_array: bool = False
    array_repeat_range: tuple[int, int] | None = None
    typical_anchor_aa: tuple[int, int] | None = None
    reference: str = ""
    notes: str = ""


@dataclass
class SystemMatch:
    """How well a locus matches one catalogue entry."""

    system_id: str
    name: str
    score: float                              # 0-1
    satisfied: list[str] = field(default_factory=list)
    unsatisfied: list[str] = field(default_factory=list)
    violations: list[str] = field(default_factory=list)

    @property
    def is_match(self) -> bool:
        return self.score >= 0.85 and not self.violations

    def explain(self) -> str:
        bits = [f"{self.name}: {self.score:.2f}"]
        bits.extend(self.shortfall())
        return "; ".join(bits)

    def shortfall(self) -> list[str]:
        """Why this locus is not that system, without restating its name."""
        bits: list[str] = []
        if self.unsatisfied:
            bits.append(f"missing {', '.join(self.unsatisfied)}")
        if self.violations:
            bits.append(f"ruled out by {', '.join(self.violations)}")
        return bits


def match_known_system(locus: Locus, system: KnownSystem) -> SystemMatch:
    """Score a locus against one described system.

    Required components drive the score; forbidden components and array
    mismatches are recorded as hard violations, which keeps a locus that merely
    resembles a system from being written off as that system.
    """
    op = predicted_operon(locus)
    labels = [g.label for g in op]
    # Neighbours outside the operon can still satisfy a requirement, but only
    # weakly, so they are checked separately.
    window_labels = [g.label for g in locus.genes]

    satisfied, unsatisfied = [], []
    for req in system.required_labels:
        if req in labels:
            satisfied.append(req)
        elif req in window_labels:
            satisfied.append(f"{req}(outside operon)")
        else:
            unsatisfied.append(req)

    violations: list[str] = []
    for bad in system.forbidden_labels:
        if bad in labels:
            violations.append(f"unexpected {bad}")

    has_array = bool(locus.arrays)
    if system.requires_array and not has_array:
        violations.append("no repeat array")
    if system.forbids_array and has_array:
        violations.append("unexpected repeat array")
    if system.array_repeat_range and has_array:
        lo, hi = system.array_repeat_range
        if not any(lo <= a.repeat_length <= hi for a in locus.arrays):
            violations.append(
                f"array repeat length {[a.repeat_length for a in locus.arrays]} "
                f"outside {lo}-{hi}"
            )
    if system.typical_anchor_aa:
        lo, hi = system.typical_anchor_aa
        n = locus.anchor.aa_length
        if n and not (lo * 0.6 <= n <= hi * 1.6):
            violations.append(f"anchor {n} aa far from typical {lo}-{hi}")

    required = len(system.required_labels) or 1
    base = len([s for s in satisfied if "(outside operon)" not in s]) / required
    weak = len([s for s in satisfied if "(outside operon)" in s]) / required
    score = min(1.0, base + 0.5 * weak)
    return SystemMatch(system.system_id, system.name, score, satisfied, unsatisfied, violations)


@dataclass
class NoveltyAssessment:
    """Whether a locus fits any described system, and how confidently."""

    best_match: SystemMatch | None
    all_matches: list[SystemMatch]
    novelty: float                 # 0 = an exact catalogue match, 1 = fits nothing
    unexplained_components: list[str] = field(default_factory=list)
    rationale: str = ""


def assess_novelty(locus: Locus, catalogue: list[KnownSystem]) -> NoveltyAssessment:
    """Compare a locus against the whole catalogue.

    Novelty is one minus the best catalogue score, but a locus that matches a
    system *and* carries a component that system never has is not fully
    explained: the unexplained component keeps some novelty on the board. That
    is exactly the situation the article's system is in -- a recognisable
    enzyme family with a component no described member has.
    """
    matches = sorted(
        (match_known_system(locus, s) for s in catalogue),
        key=lambda m: (-m.score, len(m.violations)),
    )
    if not matches:
        return NoveltyAssessment(None, [], 1.0, [], "no catalogue supplied")

    best = matches[0]
    clean = [m for m in matches if not m.violations]
    effective = max((m.score for m in clean), default=0.0)

    unexplained: list[str] = []
    if best.score >= 0.6:
        described = set(
            (best_sys := next(s for s in catalogue if s.system_id == best.system_id))
            .required_labels
        ) | set(best_sys.optional_labels)
        for g in predicted_operon(locus):
            if g.label not in described and g.label != locus.anchor.label:
                unexplained.append(f"gene:{g.label}")
        if locus.arrays and not best_sys.requires_array:
            for a in locus.arrays:
                unexplained.append(f"array:{a.copy_number}x{a.repeat_length}bp")

    novelty = 1.0 - effective
    if unexplained:
        # An unexplained component floors novelty: the locus cannot be fully
        # accounted for by the best-matching described system.
        novelty = max(novelty, 0.35 + 0.1 * min(len(unexplained), 3))
    novelty = max(0.0, min(1.0, novelty))

    if effective >= 0.85 and not unexplained:
        rationale = f"fully explained by {best.name}"
    elif effective >= 0.85:
        rationale = (
            f"matches {best.name} but carries {len(unexplained)} component(s) "
            f"that system does not have: {', '.join(unexplained)}"
        )
    elif best.score >= 0.6:
        shortfall = "; ".join(best.shortfall()) or "the match is incomplete"
        rationale = (
            f"the closest described system is {best.name}, which the gene "
            f"content matches at {best.score:.2f}, but it is {shortfall}"
        )
    else:
        rationale = (
            f"no described system scores above {best.score:.2f} "
            f"(closest: {best.name})"
        )
    return NoveltyAssessment(best, matches, novelty, unexplained, rationale)


def conservation_across_loci(loci: list[Locus]) -> dict[str, float]:
    """Fraction of loci in which each neighbour label appears.

    A partner gene that travels with the anchor across distant genomes is under
    selection to stay there. A neighbour present once is a coincidence; one
    present in most members of a clade is part of the system. The reporting
    stage quotes this number, and the triage stage treats a low value as
    grounds for demoting a "partner gene" claim.
    """
    if not loci:
        return {}
    counts: dict[str, int] = {}
    for locus in loci:
        seen = {g.label for g in predicted_operon(locus) if g.gene_id != locus.anchor.gene_id}
        if locus.arrays:
            seen.add("__array__")
        for label in seen:
            counts[label] = counts.get(label, 0) + 1
    return {label: n / len(loci) for label, n in sorted(counts.items(), key=lambda kv: -kv[1])}


def mean_intergenic_gap(loci: list[Locus]) -> float:
    """Mean anchor-to-nearest-neighbour gap across loci, in bp."""
    gaps: list[int] = []
    for locus in loci:
        ordered = sorted(locus.genes, key=lambda g: g.start)
        for i, g in enumerate(ordered):
            if g.gene_id != locus.anchor.gene_id:
                continue
            if i > 0:
                gaps.append(max(0, g.start - ordered[i - 1].end))
            if i < len(ordered) - 1:
                gaps.append(max(0, ordered[i + 1].start - g.end))
    return mean(gaps) if gaps else 0.0
