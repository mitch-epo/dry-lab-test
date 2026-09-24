"""Running a campaign end to end.

The orchestration is thin on purpose: each stage owns its own logic and writes
its own artifacts, and this module only sequences them and enforces the one
piece of cross-stage policy that matters -- the reproduce gate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from .agents import Harness, TranscriptStore, make_backend
from .campaign import CampaignConfig, Run, load_catalogue
from .rubric import Rubric
from .schemas import BenchPackage, Candidate, CandidateReport, Verdict
from .sources import open_source
from .stages import (
    run_bench,
    run_report,
    run_reproduce,
    run_search,
    run_survey,
    run_taste,
    run_triage,
)

STAGES = ["survey", "reproduce", "search", "report", "triage", "bench", "taste"]


@dataclass
class PipelineOutcome:
    run: Run
    rubric: Rubric
    gate_passed: bool = False
    stopped_at: str = ""
    candidates: list[Candidate] = field(default_factory=list)
    reports: list[CandidateReport] = field(default_factory=list)
    verdicts: list[Verdict] = field(default_factory=list)
    packages: list[BenchPackage] = field(default_factory=list)
    taste: object | None = None

    @property
    def promoted(self) -> list[Verdict]:
        return sorted(
            (v for v in self.verdicts if v.disposition == "promote"),
            key=lambda v: -v.priority,
        )


def run_campaign(
    campaign_dir: Path,
    *,
    run_id: str | None = None,
    backend: str | None = None,
    stages: list[str] | None = None,
    force: bool = False,
    labels_path: Path | None = None,
    log: Callable[[str], None] = lambda _msg: None,
) -> PipelineOutcome:
    """Execute a campaign.

    ``force`` continues past a failed reproduce gate. It exists because a
    failing gate is sometimes the thing being debugged -- but a run that used
    it is marked in the manifest, so its outputs cannot later be mistaken for
    a validated run.
    """
    campaign_dir = Path(campaign_dir)
    config = CampaignConfig.load(campaign_dir / "campaign.yaml")
    if backend:
        config.backend = backend
    wanted = stages or STAGES

    run = Run(campaign_dir / "runs" / (run_id or ""), config, run_id=run_id)
    if run_id is None:
        run = Run(campaign_dir / "runs" / run.run_id, config, run_id=run.run_id)

    rubric = run.rubric()
    source = open_source(config.resolved_source(campaign_dir))
    prov = source.provenance()

    run.manifest.source = {
        "spec": config.source,
        "kind": prov.kind,
        "description": prov.description,
        "synthetic": prov.synthetic,
        "retrieved_at": prov.retrieved_at,
        "banner": prov.banner(),
        "corpus": prov.corpus_meta,
    }
    run.manifest.synthetic_data = prov.synthetic
    run.manifest.rubric_version = rubric.version
    run.manifest.backend = config.backend
    run.save_manifest()

    transcripts = TranscriptStore(run.dir / "transcripts")
    agent_backend = make_backend(config.backend, rubric, transcripts=transcripts)
    harness = Harness(
        agent_backend, transcripts=transcripts, max_workers=config.max_workers
    )

    outcome = PipelineOutcome(run=run, rubric=rubric)
    log(f"run {run.run_id} | source {prov.kind} | backend {config.backend} | rubric v{rubric.version}")
    if prov.synthetic:
        log(prov.banner())

    # -- survey --------------------------------------------------------------
    if "survey" not in wanted:
        raise ValueError("the survey stage builds the profile every later stage needs")
    log("survey: reading literature, building the family profile and catalogue")
    survey, profile, index, systems, labeller = run_survey(run, config, source)
    log(
        f"  profile {survey.profile_length} cols "
        f"({survey.profile_core_columns} core), database {survey.database_size} proteins, "
        f"{survey.catalogue_size} described systems, "
        f"{len(labeller.models)} family profiles"
    )

    # -- reproduce (gate) ----------------------------------------------------
    if "reproduce" in wanted:
        log("reproduce: methods check against described systems")
        rep = run_reproduce(run, config, source, profile, index, systems, labeller)
        outcome.gate_passed = rep.gate_passed
        for line in rep.report().splitlines():
            if line.strip():
                log("  " + line)
        if not rep.gate_passed and not force:
            outcome.stopped_at = "reproduce"
            run.manifest.counts = {"candidates": 0, "reports": 0, "promoted": 0}
            run.save_manifest()
            log("  gate FAILED: stopping before any novel candidate is proposed")
            return outcome
        if not rep.gate_passed:
            run.manifest.stages[-1].notes.append(
                "GATE FAILED but the run continued under --force; these outputs "
                "are not validated"
            )
            run.save_manifest()
            log("  gate FAILED but continuing under --force")

    # -- search --------------------------------------------------------------
    if "search" in wanted:
        log("search: scanning for members that fit no described system")
        search, candidates = run_search(run, config, source, profile, index, systems, labeller)
        outcome.candidates = candidates
        log(f"  {search.funnel}")
        log(
            f"  {search.loci} loci -> {search.novel_loci} novel "
            f"-> {len(candidates)} candidates in {len(search.groups)} architectures"
        )
    else:
        raw = run.read_json("candidates.json") or []
        outcome.candidates = [Candidate.model_validate(c) for c in raw]

    # -- report --------------------------------------------------------------
    if "report" in wanted:
        log(f"report: writing {len(outcome.candidates)} candidate reports")
        outcome.reports = run_report(run, config, harness, rubric, outcome.candidates)
        log(f"  {len(outcome.reports)} written | {harness.stats.summary()}")
    else:
        outcome.reports = [
            CandidateReport.model_validate_json(p.read_text())
            for p in sorted((run.dir / "reports").glob("*.json"))
        ]

    # -- triage --------------------------------------------------------------
    if "triage" in wanted:
        log("triage: critically evaluating the evidence")
        outcome.verdicts = run_triage(
            run, config, harness, rubric, outcome.candidates, outcome.reports
        )
        stage = run.manifest.stage("triage")
        if stage:
            log(
                f"  promote {stage.outputs.get('promote')}, "
                f"hold {stage.outputs.get('hold')}, "
                f"eliminate {stage.outputs.get('eliminate')}"
            )
            for gate, n in sorted(
                (stage.outputs.get("gate_failures") or {}).items(), key=lambda kv: -kv[1]
            )[:6]:
                log(f"    killed by {gate}: {n}")
    else:
        raw = run.read_json("verdicts.json") or []
        outcome.verdicts = [Verdict.model_validate(v) for v in raw]

    # -- bench ---------------------------------------------------------------
    if "bench" in wanted:
        outcome.packages = run_bench(
            run, config, outcome.candidates, outcome.reports, outcome.verdicts
        )
        log(f"bench: {len(outcome.packages)} handoff packages")

    # -- taste ---------------------------------------------------------------
    if "taste" in wanted:
        log("taste: mining the report corpus for what separates promoted from set aside")
        outcome.taste = run_taste(
            run, config, harness, rubric,
            labels_path=labels_path or (campaign_dir / "labels.yaml"),
        )
        for line in outcome.taste.notes:
            log("  " + line)

    run.manifest.counts = {
        "candidates": len(outcome.candidates),
        "reports": len(outcome.reports),
        "promoted": len(outcome.promoted),
        "sessions": harness.stats.sessions,
    }
    run.save_manifest()
    run.write_text("summary.md", summarise(outcome))
    return outcome


def summarise(outcome: PipelineOutcome) -> str:
    """The run summary written to ``summary.md``."""
    run = outcome.run
    m = run.manifest
    lines = [f"# Campaign run {m.run_id}", ""]
    if m.synthetic_data:
        lines += [f"> {m.source.get('banner', '')}", ""]
    lines += [
        f"- campaign: **{m.campaign}**",
        f"- source: {m.source.get('kind')} ({m.source.get('description')})",
        f"- backend: {m.backend}",
        f"- rubric: v{m.rubric_version}",
        f"- reproduce gate: {'PASSED' if m.gate_passed else 'FAILED'}",
        "",
    ]
    if outcome.stopped_at:
        lines += [
            f"**Run stopped at the {outcome.stopped_at} stage.** The methods "
            f"check did not pass, so no novel candidate was proposed. See "
            f"`reproduce.txt`.",
            "",
        ]
        return "\n".join(lines)

    lines += ["## Stages", ""]
    for s in m.stages:
        lines.append(f"### {s.stage} ({s.duration_s:.1f}s)")
        for note in s.notes:
            lines.append(f"- {note}")
        lines.append("")

    promoted = outcome.promoted
    lines += ["## Outcome", ""]
    lines.append(
        f"{len(outcome.candidates)} candidates written up, "
        f"{len(promoted)} promoted, resolving to "
        f"{len(outcome.packages)} candidate system(s) for the bench."
    )
    lines.append("")
    if outcome.packages:
        for pkg in outcome.packages:
            lines += [
                f"### {pkg.candidate_id} -- {pkg.title}",
                "",
                f"- architecture: `{pkg.architecture}`",
                f"- independent loci: {pkg.n_instances}",
                f"- report: `reports/{pkg.candidate_id}.md`",
                f"- bench handoff: `bench/{pkg.candidate_id}.md`",
                "",
            ]
    if not promoted:
        lines += [
            "No candidate survived triage. For a family survey this is a normal "
            "result: the value of the run is the elimination, and the reasons "
            "are recorded in `verdicts.json`.",
            "",
        ]

    if outcome.taste is not None:
        lines += ["## Taste analysis", "", "```", outcome.taste.report(), "```", ""]
    return "\n".join(lines)
