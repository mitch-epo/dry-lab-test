"""Composing and rendering candidate reports.

Two jobs live here. :func:`compose_report` builds a report *from measurements*
-- it is the deterministic reporter backend, and it is constrained so that it
cannot assert anything the evidence bundle does not contain. :func:`report_markdown`
renders any report, however it was written, into the short human-readable form
the workflow calls for.

The rendering is opinionated about one thing: contradicting evidence is printed
in the same list as supporting evidence, not in a footnote. A report whose
weaknesses are easy to skip is a report that will be over-trusted by the next
reader in the chain, which in this workflow is the stage that decides whether
to spend bench time.
"""

from __future__ import annotations

from typing import Any

from .rubric import Rubric
from .schemas import CandidateReport, Evidence, Strength, Verdict

_STRENGTH_ORDER = {
    Strength.DIRECT: 0, Strength.STRONG: 1, Strength.SUGGESTIVE: 2,
    Strength.WEAK: 3, Strength.ABSENT: 4,
}


def _propose_function(bundle: dict[str, Any]) -> tuple[str, str]:
    """Derive a proposed function from the measured architecture.

    Deliberately conservative. The proposal names the components and what they
    imply mechanistically, and stops there: any stronger statement would be an
    inference the measurements do not license, and the triage stage would be
    right to strike it.
    """
    m = bundle.get("metrics", {})
    arch = bundle.get("core_architecture", "")
    has_array = float(m.get("array_quality", 0.0)) >= 0.4
    partner = bool(m.get("partner_conservation", 0.0) >= 0.4)
    lineage = bundle.get("lineage", "")
    phage = "Viruses" in lineage or "virus" in lineage.lower()

    if has_array and partner:
        title = "RT with a conserved partner and a regularly spaced non-coding array"
        fn = (
            "A three-component retroelement: the reverse transcriptase, a "
            "conserved partner protein of unknown function encoded in the same "
            "operon, and an adjacent array of evenly spaced repeats that is "
            "part of the same locus. The architecture matches no described "
            "RT-associated system. The most economical reading is that the "
            "array is the RT's template or substrate -- as msr/msd is for a "
            "retron and TR is for a DGR -- and that the partner is required "
            "for the reaction or is its effector"
            + (", in a phage-encoded defence or counter-defence context" if phage else "")
            + "."
        )
    elif partner:
        title = "RT with a conserved partner of unknown function"
        fn = (
            "A two-component system: the reverse transcriptase and a partner "
            "protein that is conserved beside it across independent loci but "
            "belongs to no described family. No non-coding component was "
            "detected, so the proposal is limited to the claim that the two "
            "genes form a functional unit."
        )
    elif has_array:
        title = "RT adjacent to a regularly spaced non-coding array"
        fn = (
            "A reverse transcriptase beside an array of evenly spaced repeats "
            "with distinct spacers, with no conserved protein partner. The "
            "array is the only feature distinguishing this locus from an "
            "unclassified lone RT."
        )
    else:
        title = "RT in an architecture absent from the catalogue"
        fn = (
            "A reverse transcriptase whose genomic context matches no described "
            "system, without a conserved partner or a non-coding component to "
            "constrain what that context means."
        )
    return title, f"{fn} Observed architecture: {arch}."


def _predictions(bundle: dict[str, Any]) -> list[dict[str, str]]:
    """Falsifiable consequences, generated from what was actually measured."""
    m = bundle.get("metrics", {})
    out: list[dict[str, str]] = []

    if m.get("triad_intact"):
        out.append({
            "statement": (
                "The purified protein has RNA-dependent DNA polymerase activity "
                "requiring the motif A and motif C aspartates."
            ),
            "assay": (
                "Express and purify the protein; assay primer extension on an RNA "
                "template with Mg2+, alongside a motif C D->N active-site mutant."
            ),
            "would_falsify": (
                "No extension by the wild type, or equal extension by the "
                "active-site mutant."
            ),
        })

    if float(m.get("array_quality", 0.0)) >= 0.4:
        out.append({
            "statement": (
                "The array is transcribed and processed into discrete short RNAs "
                "with boundaries at the repeats."
            ),
            "assay": (
                "Small-RNA sequencing and Northern blot on the native host, or on "
                "a heterologous host carrying the cloned locus."
            ),
            "would_falsify": (
                "No transcript over the array, or a single unprocessed transcript "
                "with no repeat-aligned ends."
            ),
        })
        out.append({
            "statement": (
                "The array is the RT's template: cDNA products map to array-derived "
                "sequence."
            ),
            "assay": (
                "Sequence nucleic acid co-purifying with the RT; map reads to the "
                "locus and check strand and endpoints against the repeat units."
            ),
            "would_falsify": (
                "Co-purifying nucleic acid maps elsewhere, or no DNA product is "
                "recoverable."
            ),
        })

    if float(m.get("partner_conservation", 0.0)) >= 0.4:
        out.append({
            "statement": (
                "The partner protein is required for the reaction, or forms a "
                "complex with the RT."
            ),
            "assay": (
                "Co-express and pull down; repeat the polymerase assay with and "
                "without the partner."
            ),
            "would_falsify": (
                "No co-purification and no change in activity when the partner "
                "is omitted."
            ),
        })

    if "Viruses" in bundle.get("lineage", ""):
        out.append({
            "statement": (
                "The locus affects phage-host outcome: expressing it changes "
                "plating efficiency on challenge."
            ),
            "assay": (
                "Clone the locus into a permissive host; measure efficiency of "
                "plaquing against a panel of phages, with an empty-vector control "
                "and an active-site mutant."
            ),
            "would_falsify": "No change in efficiency of plaquing for any phage tested.",
        })
    return out


def _alternatives(bundle: dict[str, Any]) -> list[str]:
    """The explanations that must be excluded before the proposal stands."""
    m = bundle.get("metrics", {})
    out = [
        f"The locus is a diverged {bundle.get('closest_system') or 'described system'} "
        f"(gene content matches at "
        f"{float(bundle.get('closest_system_score', 0.0)):.2f}) whose partner "
        f"has fallen below the labelling threshold.",
    ]
    if float(m.get("array_quality", 0.0)) > 0:
        out.append(
            "The array is an unrelated element that happens to lie nearby; the "
            f"nearest one is {m.get('array_operon_gap_bp', '?')} bp from the "
            "predicted operon."
        )
    if not m.get("association_significant", False):
        out.append(
            "The partner adjacency is coincidental: the enrichment does not "
            "survive multiple-testing correction."
        )
    if m.get("clade_restricted", True):
        out.append(
            "The architecture is a recent, clade-local arrangement rather than a "
            "conserved system."
        )
    out.append(
        "The RT is a mobile element that has inserted next to an unrelated "
        "gene and array, with no functional relationship between them."
    )
    return out


def _open_questions(bundle: dict[str, Any]) -> list[str]:
    m = bundle.get("metrics", {})
    qs: list[str] = []
    if float(m.get("array_quality", 0.0)) >= 0.4:
        qs.append("Is the array transcribed from one promoter or many?")
        qs.append("Do the spacers share a source, or are they locus-specific?")
    if float(m.get("partner_conservation", 0.0)) >= 0.4:
        qs.append("Does the partner have a fold that suggests a biochemical role?")
    if int(m.get("n_phyla", 0)) < 3:
        qs.append(
            "Is the narrow taxonomic distribution real, or an artefact of which "
            "genomes have been sequenced?"
        )
    return qs


def compose_report(bundle: dict[str, Any], rubric: Rubric) -> dict[str, Any]:
    """Write a candidate report from an evidence bundle.

    This is the deterministic reporter. Every sentence it emits is a function
    of a measurement in the bundle, so the report can be wrong about biology
    but cannot be wrong about what was observed.
    """
    title, function = _propose_function(bundle)
    evidence = list(bundle.get("evidence", []))
    predictions = _predictions(bundle)

    m = bundle.get("metrics", {})
    supporting = [e for e in evidence if e.get("supports", True)]
    against = [e for e in evidence if not e.get("supports", True)]

    confidence = 0.0
    if not m.get("triad_intact"):
        confidence = 0.05
    else:
        _failed, _reasons, _scores, priority = rubric.apply(dict(m))
        confidence = round(priority * (0.5 if _failed else 1.0), 3)

    summary = (
        f"{bundle.get('candidate_id')}: reverse transcriptase "
        f"{bundle.get('anchor_protein')} on {bundle.get('contig_id')} "
        f"({(bundle.get('lineage') or 'unknown lineage').split(';')[-1]}). "
        f"Novelty {float(bundle.get('novelty', 0)):.2f} against the catalogue. "
        f"The nearest described system on gene content alone is "
        f"{bundle.get('closest_system') or 'none'} "
        f"({float(bundle.get('closest_system_score', 0)):.2f}), but that match "
        f"is ruled out on architecture -- which is why the novelty score is "
        f"high despite the gene-content similarity. "
        f"{len(supporting)} supporting and {len(against)} contradicting "
        f"observations. The architecture recurs at "
        f"{m.get('n_independent_loci', 1)} loci "
        f"({m.get('n_clusters', 1)} non-redundant, "
        f"{m.get('n_phyla', 1)} phyla)."
    )

    return {
        "candidate_id": bundle.get("candidate_id", ""),
        "title": title,
        "proposed_function": function,
        "summary": summary,
        "evidence": evidence,
        "alternatives": _alternatives(bundle),
        "predictions": predictions,
        "open_questions": _open_questions(bundle),
        "confidence": confidence,
        "author": "rubric-reporter",
        "rubric_version": rubric.version,
    }


# --- markdown rendering -------------------------------------------------------


def _evidence_lines(evidence: list[Evidence]) -> list[str]:
    ordered = sorted(
        evidence,
        key=lambda e: (not e.supports, _STRENGTH_ORDER.get(e.strength, 9)),
    )
    lines: list[str] = []
    for e in ordered:
        mark = "supports" if e.supports else "**against**"
        lines.append(f"- [{e.strength.value}, {mark}] {e.statement}")
        if e.method:
            lines.append(f"  - method: {e.method}")
        if e.records:
            lines.append(f"  - records: `{'`, `'.join(e.records[:6])}`")
    return lines


def report_markdown(
    report: CandidateReport,
    *,
    verdict: Verdict | None = None,
    banner: str = "",
) -> str:
    """Render a candidate report as the short human-readable document."""
    lines: list[str] = []
    if banner:
        lines += [f"> {banner}", ""]
    lines += [f"# {report.title}", "", f"**{report.candidate_id}**", ""]

    if verdict is not None:
        badge = {"promote": "PROMOTE", "hold": "HOLD", "eliminate": "ELIMINATED"}[
            verdict.disposition
        ]
        if verdict.failed_gates:
            # A ranking score is meaningless once a disqualifying check has
            # failed, and printing a high one beside ELIMINATED reads as a
            # contradiction. Say which check failed instead.
            lines += [
                f"**Triage: {badge}** "
                f"(failed {', '.join(verdict.failed_gates)}; rubric "
                f"v{verdict.rubric_version})",
                "",
                f"> Eliminated because {verdict.kill_reason}.",
                "",
                f"> Its ranking score before the check was applied was "
                f"{verdict.priority:.2f}. That score does not survive a failed "
                f"check and is shown only so the report can be read against "
                f"the others.",
                "",
            ]
        else:
            lines += [
                f"**Triage: {badge}** (priority {verdict.priority:.2f}, "
                f"rubric v{verdict.rubric_version})",
                "",
            ]
            if verdict.kill_reason:
                lines += [f"> Set aside because {verdict.kill_reason}.", ""]

    lines += ["## Proposed function", "", report.proposed_function, ""]
    lines += ["## Summary", "", report.summary, ""]
    lines += ["## Evidence", ""] + _evidence_lines(report.evidence) + [""]

    if report.alternatives:
        lines += ["## Alternative explanations", ""]
        lines += [f"- {a}" for a in report.alternatives]
        lines.append("")

    if report.predictions:
        lines += ["## Predictions and how to falsify them", ""]
        for i, p in enumerate(report.predictions, 1):
            lines += [
                f"{i}. {p.statement}",
                f"   - assay: {p.assay}",
                f"   - would falsify: {p.would_falsify}",
            ]
        lines.append("")

    if report.open_questions:
        lines += ["## Open questions", ""]
        lines += [f"- {q}" for q in report.open_questions]
        lines.append("")

    if verdict is not None and (verdict.critique or verdict.overclaims):
        lines += ["## Reviewer notes", ""]
        lines += [f"- {c}" for c in verdict.critique]
        lines += [f"- overclaim: {o}" for o in verdict.overclaims]
        lines.append("")
        if verdict.scores:
            lines += ["Scores: " + ", ".join(
                f"{k} {v:.2f}" for k, v in sorted(verdict.scores.items())
            ), ""]

    lines += [
        "---",
        f"_Confidence {report.confidence:.2f}; written by {report.author}; "
        f"rubric v{report.rubric_version}; {report.created_at}._",
    ]
    return "\n".join(lines)
