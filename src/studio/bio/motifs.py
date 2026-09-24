"""Motif models: PROSITE-style patterns with spacing constraints.

A profile hit says a protein belongs to a family. It does not say the protein
can still do the family's chemistry -- a pseudogene, a fragment or a
catalytically dead paralogue all score well against the profile. The workflow's
first real filter is whether the catalytic residues are present and correctly
spaced, so motifs are modelled explicitly rather than left implicit in the PSSM.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from .seq import clean_protein

_RANGE = re.compile(r"^x\((\d+)(?:,(\d+))?\)$", re.IGNORECASE)


def prosite_to_regex(pattern: str) -> str:
    """Compile a PROSITE-style pattern into a Python regular expression.

    Supported syntax: ``[ABC]`` any-of, ``{ABC}`` none-of, ``x`` any residue,
    ``x(n)`` / ``x(n,m)`` runs, elements separated by ``-``, and a trailing
    ``>`` / leading ``<`` for terminal anchors.
    """
    out: list[str] = []
    pattern = pattern.strip()
    if pattern.startswith("<"):
        out.append("^")
        pattern = pattern[1:]
    anchor_end = pattern.endswith(">")
    if anchor_end:
        pattern = pattern[:-1]

    for element in (e for e in pattern.split("-") if e):
        m = _RANGE.match(element)
        if m:
            lo, hi = m.group(1), m.group(2)
            out.append(f"[A-Z]{{{lo},{hi}}}" if hi else f"[A-Z]{{{lo}}}")
        elif element.startswith("[") and element.endswith("]"):
            out.append(f"[{element[1:-1].upper()}]")
        elif element.startswith("{") and element.endswith("}"):
            out.append(f"[^{element[1:-1].upper()}]")
        elif element.lower() == "x":
            out.append("[A-Z]")
        else:
            out.append(re.escape(element.upper()))
    if anchor_end:
        out.append("$")
    return "".join(out)


def pattern_information(pattern: str) -> float:
    """Bits of information in a PROSITE-style pattern.

    Each fixed or restricted position contributes ``log2(20 / |allowed set|)``;
    wildcards contribute nothing. This measures how surprising a match is, and
    the scanner uses it to prefer an informative placement over a permissive
    one -- a ``Y-x-D-D`` match is far rarer than ``L-x-D-D`` and should win when
    both are geometrically acceptable.
    """
    bits = 0.0
    pattern = pattern.strip().lstrip("<").rstrip(">")
    for element in (e for e in pattern.split("-") if e):
        if _RANGE.match(element) or element.lower() == "x":
            continue
        if element.startswith("[") and element.endswith("]"):
            allowed = len(set(element[1:-1].upper())) or 20
        elif element.startswith("{") and element.endswith("}"):
            allowed = max(1, 20 - len(set(element[1:-1].upper())))
        else:
            allowed = 1
            bits += len(element) * math.log2(20.0)
            continue
        bits += math.log2(20.0 / allowed)
    return bits


@dataclass
class MotifSpec:
    """A named motif: a pattern plus what it means if it is missing."""

    name: str
    pattern: str
    description: str = ""
    essential: bool = False       # absence is disqualifying, not merely notable
    min_offset: int | None = None  # earliest plausible start, in residues
    max_offset: int | None = None

    def __post_init__(self) -> None:
        self._rx = re.compile(prosite_to_regex(self.pattern))
        self.information = pattern_information(self.pattern)

    def find(self, seq: str) -> list[tuple[int, str]]:
        """All matches as ``(start, matched_text)``, including overlaps."""
        seq = clean_protein(seq)
        out: list[tuple[int, str]] = []
        pos = 0
        while pos < len(seq):
            m = self._rx.search(seq, pos)
            if not m:
                break
            out.append((m.start(), m.group(0)))
            pos = m.start() + 1
        return out


@dataclass
class MotifHit:
    name: str
    start: int
    text: str
    essential: bool
    information: float = 0.0


@dataclass
class MotifReport:
    """Outcome of scanning one protein with a :class:`MotifSet`."""

    hits: list[MotifHit] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    missing_essential: list[str] = field(default_factory=list)
    spacing_violations: list[str] = field(default_factory=list)
    order_ok: bool = True

    @property
    def intact(self) -> bool:
        """True when every essential motif is present, ordered and correctly spaced."""
        return (
            not self.missing_essential
            and not self.spacing_violations
            and self.order_ok
        )

    @property
    def completeness(self) -> float:
        total = len(self.hits) + len(self.missing)
        return len(self.hits) / total if total else 0.0

    def summary(self) -> str:
        found = ", ".join(f"{h.name}@{h.start}:{h.text}" for h in self.hits) or "none"
        parts = [f"found[{found}]"]
        if self.missing:
            parts.append(f"missing[{','.join(self.missing)}]")
        if self.spacing_violations:
            parts.append(f"spacing[{'; '.join(self.spacing_violations)}]")
        if not self.order_ok:
            parts.append("order[out of order]")
        return " ".join(parts)


@dataclass
class MotifSet:
    """An ordered set of motifs with allowed inter-motif spacing.

    ``spacing`` maps ``(earlier, later)`` motif names to an inclusive
    ``(min, max)`` residue gap between the end of the first and the start of
    the second.
    """

    name: str
    motifs: list[MotifSpec]
    spacing: dict[tuple[str, str], tuple[int, int]] = field(default_factory=dict)
    ordered: bool = True

    max_occurrences: int = 40

    def _gap_ok(self, a: MotifHit, b: MotifHit) -> bool:
        """Check any spacing constraint registered for the pair, in either order."""
        for (first, second), (lo, hi) in self.spacing.items():
            if a.name == first and b.name == second:
                gap = b.start - (a.start + len(a.text))
            elif a.name == second and b.name == first:
                gap = a.start - (b.start + len(b.text))
            else:
                continue
            if gap < lo or gap > hi:
                return False
        return True

    def scan(self, seq: str) -> MotifReport:
        """Find the best mutually consistent placement of the whole motif set.

        Choosing each motif's first match independently is unsafe: a short,
        permissive pattern matching by chance early in the sequence can block
        the genuine downstream match of a later motif, and the set then reads as
        disrupted when it is intact. Instead every occurrence of every motif is
        enumerated and the assignment satisfying the ordering and spacing
        constraints is searched for directly, with each motif also allowed to be
        absent. Assignments are ranked by essential motifs placed, then total
        motifs placed, then compactness -- a real catalytic core is contiguous.
        """
        seq = clean_protein(seq)

        options: list[list[MotifHit | None]] = []
        for spec in self.motifs:
            occ: list[MotifHit | None] = []
            for start, text in spec.find(seq)[: self.max_occurrences]:
                if spec.min_offset is not None and start < spec.min_offset:
                    continue
                if spec.max_offset is not None and start > spec.max_offset:
                    continue
                occ.append(MotifHit(spec.name, start, text, spec.essential, spec.information))
            occ.append(None)  # absence is always an option
            options.append(occ)

        best: list[MotifHit | None] | None = None
        best_key: tuple[int, int, float, int] | None = None
        remaining_essential = [
            sum(1 for s in self.motifs[i:] if s.essential) for i in range(len(self.motifs) + 1)
        ]

        def key_of(assign: list[MotifHit | None]) -> tuple[int, int, float, int]:
            hits = [h for h in assign if h is not None]
            ess = sum(1 for h in hits if h.essential)
            info = sum(h.information for h in hits)
            span = (max(h.start for h in hits) - min(h.start for h in hits)) if hits else 0
            return (ess, len(hits), info, -span)

        def search(i: int, chosen: list[MotifHit | None], placed: list[MotifHit]) -> None:
            nonlocal best, best_key
            if i == len(self.motifs):
                k = key_of(chosen)
                if best_key is None or k > best_key:
                    best, best_key = list(chosen), k
                return
            # Prune: even placing every remaining essential motif cannot beat the best.
            if best_key is not None:
                ceiling = sum(1 for h in placed if h.essential) + remaining_essential[i]
                if ceiling < best_key[0]:
                    return
            for cand in options[i]:
                if cand is not None:
                    if self.ordered and placed and cand.start < placed[-1].start:
                        continue
                    if any(not self._gap_ok(p, cand) for p in placed):
                        continue
                    chosen.append(cand)
                    placed.append(cand)
                    search(i + 1, chosen, placed)
                    chosen.pop()
                    placed.pop()
                else:
                    chosen.append(None)
                    search(i + 1, chosen, placed)
                    chosen.pop()

        search(0, [], [])
        assignment = best or [None] * len(self.motifs)

        report = MotifReport()
        placed_by_name: dict[str, MotifHit] = {}
        for spec, hit in zip(self.motifs, assignment):
            if hit is None:
                report.missing.append(spec.name)
                if spec.essential:
                    report.missing_essential.append(spec.name)
            else:
                report.hits.append(hit)
                placed_by_name[spec.name] = hit

        # Spacing constraints are satisfied by construction for placed pairs; a
        # constraint is reported as violated only when both motifs were placed
        # and the search had to break it to place them at all (impossible here),
        # so what remains worth reporting is a constraint left unevaluated.
        for (first, second), (lo, hi) in self.spacing.items():
            a, b = placed_by_name.get(first), placed_by_name.get(second)
            if a is None or b is None:
                continue
            gap = b.start - (a.start + len(a.text))
            if gap < lo or gap > hi:
                report.spacing_violations.append(
                    f"{first}->{second} gap {gap} outside [{lo},{hi}]"
                )

        if self.ordered:
            order = [h.start for h in report.hits]
            report.order_ok = order == sorted(order)
        return report


# --- Reverse transcriptase catalytic core ------------------------------------
#
# The RT palm domain carries a two-metal-ion active site formed by an aspartate
# in motif A and the aspartate pair of the motif C YxDD signature. Motif B
# contributes the conserved basic residue that positions the incoming dNTP.
# These are the classical Xiong-Eickbush RT motifs; the spacing bounds are
# generous because the studio must not exclude a genuinely divergent family
# member on spacing alone.

RT_MOTIFS = MotifSet(
    name="RT-palm",
    motifs=[
        MotifSpec(
            name="A",
            pattern="[LIVMF]-x-[DN]-[LIVMFYAG]",
            description="Motif A: catalytic aspartate in a hydrophobic bed",
            essential=True,
        ),
        MotifSpec(
            name="B",
            pattern="[KR]-x(2,4)-[LIVMFYA]-x-[LIVMFYAG]-[KRQ]",
            description="Motif B: dNTP-positioning basic residue",
            essential=False,
        ),
        MotifSpec(
            name="C",
            pattern="[YF]-x-[DN]-[DN]",
            description="Motif C: YxDD catalytic aspartate pair",
            essential=True,
        ),
        MotifSpec(
            name="E",
            pattern="[LIVMFY]-[GAST]-x-[LIVMFY]-x-[KRQNH]",
            description="Motif E: primer grip",
            essential=False,
        ),
    ],
    spacing={
        # Motif A sits well upstream of the YxDD pair in every characterised RT
        # clade. The window is wide because divergent families insert loops
        # between them, but it is tight enough to reject the chance
        # hydrophobic-x-D matches that a four-residue pattern inevitably picks
        # up: motif A on its own is never treated as evidence, only motif A at
        # a plausible distance from a real motif C.
        ("A", "C"): (40, 180),
        ("C", "E"): (0, 90),
    },
)


def catalytic_triad_intact(seq: str) -> tuple[bool, str]:
    """Check the RT two-metal-ion triad: motif A aspartate plus the motif C pair.

    This asks its own question rather than reading off :meth:`MotifSet.scan`.
    The general scan optimises the placement of the whole set, and the
    assignment that places the most motifs is not always the one containing the
    genuine catalytic pair -- a chance ``YxDN`` elsewhere in the protein can win
    that competition and make an intact enzyme read as dead. Here every motif A
    and motif C occurrence is considered, and the protein counts as intact if
    *any* correctly spaced pair carries both aspartates.

    Returns ``(intact, explanation)``. The explanation names residues and
    positions so it can be quoted verbatim in a candidate report.
    """
    seq = clean_protein(seq)
    spec_a = next(s for s in RT_MOTIFS.motifs if s.name == "A")
    spec_c = next(s for s in RT_MOTIFS.motifs if s.name == "C")
    lo, hi = RT_MOTIFS.spacing[("A", "C")]

    a_hits = spec_a.find(seq)
    c_hits = spec_c.find(seq)

    if not c_hits:
        return False, "no motif C (YxDD-like) match; catalytic aspartate pair absent"

    spaced: list[tuple[tuple[int, str], tuple[int, str]]] = []
    for a_start, a_text in a_hits:
        for c_start, c_text in c_hits:
            gap = c_start - (a_start + len(a_text))
            if lo <= gap <= hi:
                spaced.append(((a_start, a_text), (c_start, c_text)))

    if not spaced:
        best_c = c_hits[0]
        if not a_hits:
            return False, (
                f"motif C present ({best_c[1]}@{best_c[0] + 1}) but no motif A "
                f"aspartate found upstream"
            )
        return False, (
            f"motif A and motif C both present but never at a plausible "
            f"separation ({lo}-{hi} residues): closest pair is "
            f"{a_hits[0][1]}@{a_hits[0][0] + 1} / {best_c[1]}@{best_c[0] + 1}"
        )

    for (a_start, a_text), (c_start, c_text) in spaced:
        d_pos = [c_start + i for i, ch in enumerate(c_text) if ch == "D"]
        a_d = [a_start + i for i, ch in enumerate(a_text) if ch == "D"]
        if len(d_pos) >= 2 and a_d:
            return True, (
                f"catalytic triad intact: D{a_d[0] + 1} (motif A {a_text}"
                f"@{a_start + 1}), D{d_pos[0] + 1}/D{d_pos[1] + 1} "
                f"(motif C {c_text}@{c_start + 1})"
            )

    # Correctly spaced pairs exist but none carries a complete triad.
    (a_start, a_text), (c_start, c_text) = spaced[0]
    d_pos = [i for i, ch in enumerate(c_text) if ch == "D"]
    if len(d_pos) < 2:
        return False, (
            f"motif C at {c_start + 1} reads {c_text}: the aspartate pair is "
            f"substituted (N for D), consistent with a catalytically dead copy"
        )
    return False, (
        f"motif A at {a_start + 1} reads {a_text}: the position that should "
        f"carry the catalytic aspartate is substituted"
    )
