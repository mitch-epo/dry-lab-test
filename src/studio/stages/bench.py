"""Stage 6 -- bench: the wet-lab handoff for candidates that survived.

"When a candidate survives our review, we test it in the laboratory, expressing
the protein in standard laboratory strains and characterizing it biochemically
and structurally."

This stage produces the document that crosses from dry to wet: what to build,
what to measure, which controls make the measurement interpretable, and --
stated first, because it is the part that gets dropped -- what result would
kill the hypothesis.

The package is deliberately mundane. Expression constructs, an active-site
mutant, an empty-vector control, a defined readout. Nothing here designs an
organism or proposes a novel agent; it is the standard characterisation any
lab would run on a protein of unknown function, and the biosafety note records
the scope the handoff assumes.
"""

from __future__ import annotations

import time
from pathlib import Path

from ..campaign import CampaignConfig, Run
from ..schemas import BenchPackage, Candidate, CandidateReport, StageRecord, Verdict


def _constructs(candidate: Candidate, report: CandidateReport) -> list[dict[str, str]]:
    m = candidate.metrics
    out = [
        {
            "name": "pEXP-anchor-WT",
            "insert": f"{candidate.anchor_protein} coding sequence, codon-optimised",
            "tag": "N-terminal His6-SUMO, cleavable",
            "host": "E. coli BL21(DE3)",
            "purpose": "purify the enzyme for in vitro polymerase assays",
        },
        {
            "name": "pEXP-anchor-DD/NN",
            "insert": "same, with the motif C aspartate pair mutated to asparagine",
            "tag": "N-terminal His6-SUMO, cleavable",
            "host": "E. coli BL21(DE3)",
            "purpose": (
                "active-site control: every activity attributed to the enzyme "
                "must disappear here, or it is not the enzyme's activity"
            ),
        },
    ]
    if float(m.get("partner_conservation", 0.0)) >= 0.4:
        out.append({
            "name": "pEXP-partner",
            "insert": "the conserved partner gene from the same locus",
            "tag": "C-terminal Strep-II",
            "host": "E. coli BL21(DE3)",
            "purpose": "co-expression and pull-down; test whether it is required for activity",
        })
    if float(m.get("array_quality", 0.0)) >= 0.4:
        out.append({
            "name": "pLOC-full",
            "insert": (
                f"the whole locus including the {m.get('array_copies', 0)}-copy "
                f"array ({m.get('array_repeat_len', 0)} bp repeats), native promoter"
            ),
            "tag": "none",
            "host": "E. coli MG1655 and a permissive native-adjacent host",
            "purpose": "test array transcription and processing in context",
        })
        out.append({
            "name": "pLOC-array-deleted",
            "insert": "the same locus with the array removed",
            "tag": "none",
            "host": "as above",
            "purpose": (
                "the array's contribution is only interpretable against its absence"
            ),
        })
    return out


def _assays(candidate: Candidate, report: CandidateReport) -> list[dict[str, str]]:
    m = candidate.metrics
    out = [
        {
            "name": "RNA-dependent DNA polymerase assay",
            "readout": "primer extension on a labelled RNA template, denaturing PAGE",
            "conditions": "Mg2+ and Mn2+ titration, 30 and 37 C, +/- RNase H",
            "decides": "whether the protein is an active RT at all",
        },
        {
            "name": "Size-exclusion chromatography with multi-angle light scattering",
            "readout": "oligomeric state of the enzyme alone and with the partner",
            "conditions": "physiological salt",
            "decides": "whether enzyme and partner form a defined complex",
        },
    ]
    if float(m.get("array_quality", 0.0)) >= 0.4:
        out.append({
            "name": "Small-RNA sequencing and Northern blot",
            "readout": "size and boundaries of transcripts from the array",
            "conditions": "native host if available, otherwise the cloned locus",
            "decides": (
                "whether the array is expressed as discrete short RNAs, as the "
                "proposal requires"
            ),
        })
        out.append({
            "name": "Co-purifying nucleic acid sequencing",
            "readout": "identity, strand and endpoints of nucleic acid bound to the enzyme",
            "conditions": "pull-down from the host carrying the full locus",
            "decides": "whether the array is the template the enzyme copies",
        })
    if "Viruses" in candidate.lineage:
        out.append({
            "name": "Efficiency of plaquing",
            "readout": "plaque counts on hosts carrying the locus versus empty vector",
            "conditions": "panel of coliphages, 30 and 37 C",
            "decides": "whether the locus has a phage-defence phenotype",
        })
    out.append({
            "name": "Structure determination",
            "readout": "cryo-EM or crystallography of the enzyme, and of the complex",
            "conditions": "with and without the partner; with nucleic acid where it co-purifies",
            "decides": "how the partner and the nucleic acid are positioned relative to the active site",
    })
    return out


def _falsification(candidate: Candidate, report: CandidateReport) -> list[str]:
    """What would kill the hypothesis. Listed first in the rendered handoff."""
    out = [
        "The purified protein shows no primer extension on any template under "
        "any tested condition: the family assignment is wrong or the protein "
        "is not an active enzyme.",
        "The active-site mutant retains the activity: the readout is a "
        "contaminant, not the enzyme.",
    ]
    m = candidate.metrics
    if float(m.get("array_quality", 0.0)) >= 0.4:
        out.append(
            "No transcript is detectable over the array, or the transcript is "
            "a single unprocessed species: the array is not the non-coding "
            "component the proposal requires."
        )
        out.append(
            "Nucleic acid co-purifying with the enzyme does not map to the "
            "array: the array is adjacent but not the template."
        )
    if float(m.get("partner_conservation", 0.0)) >= 0.4:
        out.append(
            "The partner neither co-purifies nor changes activity: the "
            "conserved adjacency does not reflect a physical or functional "
            "partnership."
        )
    out.extend(p.would_falsify for p in report.predictions)
    # De-duplicate while preserving order.
    seen: set[str] = set()
    return [x for x in out if not (x in seen or seen.add(x))]


def build_package(
    candidate: Candidate,
    report: CandidateReport,
    *,
    supporting: list[str] | None = None,
) -> BenchPackage:
    supporting = supporting or []
    return BenchPackage(
        candidate_id=candidate.candidate_id,
        title=report.title,
        proposed_function=report.proposed_function,
        architecture=candidate.core_architecture,
        supporting_loci=supporting,
        n_instances=len(supporting) + 1,
        constructs=_constructs(candidate, report),
        assays=_assays(candidate, report),
        controls=[
            "empty vector in the same host and induction conditions",
            "catalytically dead active-site mutant for every activity readout",
            "a described RT from the catalogue (a retron RT) as a positive control "
            "for the polymerase assay",
            "mock pull-down from a host not expressing the tagged protein",
        ],
        falsification=_falsification(candidate, report),
        biosafety_notes=[
            "Standard BSL-1 work: heterologous expression of a protein of unknown "
            "function in a laboratory E. coli strain, and characterisation in vitro.",
            "Phage-challenge experiments use established laboratory coliphages and "
            "laboratory host strains.",
            "The handoff covers characterisation only. Any result suggesting a "
            "phenotype beyond the scope above should be reviewed before follow-up.",
        ],
    )


def package_markdown(pkg: BenchPackage, *, banner: str = "") -> str:
    lines: list[str] = []
    if banner:
        lines += [f"> {banner}", ""]
    lines += [
        f"# Bench handoff -- {pkg.candidate_id}",
        "",
        f"**{pkg.title}**",
        "",
        f"Architecture `{pkg.architecture}`, seen at {pkg.n_instances} "
        f"independent {'locus' if pkg.n_instances == 1 else 'loci'}."
        + (
            " This handoff covers the system; the other instances are listed "
            "at the end and are not separate experiments."
            if pkg.supporting_loci else ""
        ),
        "",
        "## Proposed function",
        "",
        pkg.proposed_function,
        "",
        "## What would kill this hypothesis",
        "",
    ]
    lines += [f"- {f}" for f in pkg.falsification]
    lines += ["", "## Constructs", ""]
    for c in pkg.constructs:
        lines.append(f"### {c['name']}")
        lines += [f"- {k}: {v}" for k, v in c.items() if k != "name"]
        lines.append("")
    lines += ["## Assays", ""]
    for a in pkg.assays:
        lines.append(f"### {a['name']}")
        lines += [f"- {k}: {v}" for k, v in a.items() if k != "name"]
        lines.append("")
    lines += ["## Controls", ""] + [f"- {c}" for c in pkg.controls]
    lines += ["", "## Scope and biosafety", ""] + [f"- {n}" for n in pkg.biosafety_notes]
    if pkg.supporting_loci:
        lines += [
            "",
            "## Other loci with this architecture",
            "",
            "Independent instances supporting the same proposal. They are "
            "evidence for the system, not additional experiments.",
            "",
        ]
        lines += [f"- `{c}`" for c in pkg.supporting_loci]
    lines.append("")
    return "\n".join(lines)


def run_bench(
    run: Run,
    config: CampaignConfig,
    candidates: list[Candidate],
    reports: list[CandidateReport],
    verdicts: list[Verdict],
) -> list[BenchPackage]:
    """Build a handoff for every promoted candidate."""
    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    by_cand = {c.candidate_id: c for c in candidates}
    by_report = {r.candidate_id: r for r in reports}
    promoted = [v for v in verdicts if v.disposition == "promote"]
    promoted.sort(key=lambda v: -v.priority)

    # One handoff per architecture, not per locus. Fifteen loci sharing an
    # architecture are fifteen instances of one candidate system, and the
    # decision the bench is being asked to make is about the system. Reporting
    # them separately would inflate the apparent yield of the survey and send
    # the same experiment to the bench fifteen times.
    grouped: dict[str, list[Verdict]] = {}
    for v in promoted:
        c = by_cand.get(v.candidate_id)
        if c is None:
            continue
        grouped.setdefault(c.group_signature or c.core_architecture, []).append(v)

    bench_dir = Path(run.dir) / "bench"
    banner = run.manifest.source.get("banner", "")
    packages: list[BenchPackage] = []
    for signature, members in sorted(
        grouped.items(), key=lambda kv: -max(v.priority for v in kv[1])
    ):
        lead = members[0]
        c, r = by_cand.get(lead.candidate_id), by_report.get(lead.candidate_id)
        if c is None or r is None:
            continue
        pkg = build_package(
            c, r, supporting=[m.candidate_id for m in members[1:]]
        )
        packages.append(pkg)
        (bench_dir / f"{pkg.candidate_id}.json").write_text(pkg.model_dump_json(indent=1))
        (bench_dir / f"{pkg.candidate_id}.md").write_text(
            package_markdown(pkg, banner=banner)
        )

    notes = [
        f"{len(packages)} candidate system(s) from {len(promoted)} promoted loci "
        f"({len(grouped)} distinct architectures)"
    ]
    if not packages:
        notes.append(
            "nothing was promoted, so there is nothing to hand to the bench: "
            "the survey ends here"
        )

    run.record_stage(StageRecord(
        stage="bench",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"promoted": len(promoted)},
        outputs={"packages": len(packages), "architectures": len(grouped)},
        notes=notes,
    ))
    return packages
