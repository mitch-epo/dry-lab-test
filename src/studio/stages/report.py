"""Stage 4 -- report: one short, human-readable write-up per candidate.

"...and writes a short, human-readable report for each candidate that proposes
a function and describes the evidence supporting its claims."

One session per candidate, run in parallel by the harness. Each session
receives only that candidate's evidence bundle, so a report cannot borrow
confidence from a sibling candidate it never saw -- which is also what makes
the reports independent enough for the triage and taste stages to treat as
separate observations.

The report stage does not decide anything. It proposes and it cites. Deciding
is the next stage's job, and keeping the two apart is what lets the reviewer
disagree with a report on the evidence rather than on the prose.
"""

from __future__ import annotations

import time
from pathlib import Path

from ..agents import Harness, Session
from ..campaign import CampaignConfig, Run
from ..render import report_markdown
from ..rubric import Rubric
from ..schemas import Candidate, CandidateReport, Evidence, Prediction, StageRecord

_SCHEMA_HINT = """{
  "candidate_id": "string",
  "title": "one line naming the architecture, not a conclusion",
  "proposed_function": "what the locus does and why the architecture implies it",
  "summary": "three or four sentences a reviewer can read in isolation",
  "evidence": [{"key": "...", "statement": "...", "value": null,
                "strength": "direct|strong|suggestive|weak|absent",
                "method": "...", "records": ["..."], "supports": true}],
  "alternatives": ["explanations that must be excluded before the proposal stands"],
  "predictions": [{"statement": "...", "assay": "...", "would_falsify": "..."}],
  "open_questions": ["..."],
  "confidence": 0.0
}"""

_INSTRUCTIONS = """You are writing a candidate report in a protein-family survey.

Your job is to propose a function for one locus and to describe the evidence
for it accurately. You are not deciding whether the candidate is worth testing;
a separate reviewer does that, and will check every claim you make against the
measurements in the bundle.

Rules:
- Every claim must cite a measurement from the bundle. If a number is not in
  the bundle, you do not have it.
- Carry contradicting evidence into the report with the same prominence as
  supporting evidence. Do not bury it.
- Assign each piece of evidence a strength that the measurement actually
  supports. An association with p = 0.04 at four loci is suggestive, not strong.
- State the alternative explanations, and for each, what observation would
  distinguish it from your proposal.
- Give predictions that one experiment could falsify. A prediction that cannot
  fail is not a prediction.
- Prefer "the architecture matches no described system" over "this is a novel
  system". The first is what was measured.

{rubric}
"""


def _bundle(candidate: Candidate) -> dict:
    """The evidence a reporter session sees. Nothing else is in scope."""
    return {
        "candidate_id": candidate.candidate_id,
        "anchor_protein": candidate.anchor_protein,
        "contig_id": candidate.contig_id,
        "lineage": candidate.lineage,
        "architecture": candidate.architecture,
        "core_architecture": candidate.core_architecture,
        "novelty": candidate.novelty,
        "closest_system": candidate.closest_system,
        "closest_system_score": candidate.closest_system_score,
        "metrics": candidate.metrics,
        "evidence": [e.model_dump(mode="json") for e in candidate.evidence],
        "homolog_loci": candidate.homolog_loci,
    }


def _coerce(payload: dict, candidate: Candidate, rubric: Rubric) -> CandidateReport:
    """Validate a backend's output into a report, tolerating a loose model reply."""
    evidence: list[Evidence] = []
    for e in payload.get("evidence", []):
        if isinstance(e, Evidence):
            evidence.append(e)
            continue
        try:
            evidence.append(Evidence.model_validate(e))
        except Exception:
            # A malformed evidence item is dropped rather than guessed at, and
            # the loss is visible because the report's evidence count changes.
            continue
    if not evidence:
        evidence = list(candidate.evidence)

    predictions: list[Prediction] = []
    for p in payload.get("predictions", []):
        try:
            predictions.append(
                p if isinstance(p, Prediction) else Prediction.model_validate(p)
            )
        except Exception:
            continue

    return CandidateReport(
        candidate_id=payload.get("candidate_id") or candidate.candidate_id,
        title=payload.get("title") or candidate.core_architecture,
        proposed_function=payload.get("proposed_function", ""),
        summary=payload.get("summary", ""),
        evidence=evidence,
        alternatives=list(payload.get("alternatives", [])),
        predictions=predictions,
        open_questions=list(payload.get("open_questions", [])),
        confidence=float(payload.get("confidence", 0.0) or 0.0),
        author=payload.get("author", ""),
        rubric_version=int(payload.get("rubric_version", rubric.version)),
    )


def run_report(
    run: Run,
    config: CampaignConfig,
    harness: Harness,
    rubric: Rubric,
    candidates: list[Candidate],
) -> list[CandidateReport]:
    """Write one report per candidate, in parallel."""
    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    instructions = _INSTRUCTIONS.format(rubric=rubric.as_prompt())
    sessions = [
        Session(
            session_id=f"report.{c.candidate_id}",
            role="reporter",
            instructions=instructions,
            task=(
                "Write the candidate report for the locus in this bundle. "
                "Propose a function, describe the evidence for and against it, "
                "list the alternatives, and give falsifiable predictions."
            ),
            bundle=_bundle(c),
            schema_hint=_SCHEMA_HINT,
            metadata={"candidate_id": c.candidate_id},
        )
        for c in candidates
    ]

    results = harness.run(sessions)
    by_id = {c.candidate_id: c for c in candidates}

    reports: list[CandidateReport] = []
    failed: list[str] = []
    for session, res in zip(sessions, results):
        cid = session.metadata["candidate_id"]
        if not res.ok:
            failed.append(f"{cid}: {res.error}")
            continue
        report = _coerce(res.output, by_id[cid], rubric)
        if not report.author:
            report.author = res.backend
        reports.append(report)
        # The number of falsifiable predictions feeds the testability score, so
        # it has to travel back onto the candidate before triage reads it.
        by_id[cid].metrics["n_predictions"] = len(report.predictions)

    reports_dir = Path(run.dir) / "reports"
    for r in reports:
        (reports_dir / f"{r.candidate_id}.json").write_text(r.model_dump_json(indent=1))

    run.write_json("candidates.json", [c.model_dump() for c in candidates])
    run.record_stage(StageRecord(
        stage="report",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"candidates": len(candidates)},
        outputs={"reports": len(reports), "failed": len(failed)},
        sessions=len(sessions),
        notes=[harness.stats.summary()] + ([f"failed: {'; '.join(failed[:5])}"] if failed else []),
    ))
    return reports
