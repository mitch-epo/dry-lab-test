"""Stage 5 -- triage: critically evaluate the evidence, and eliminate.

"In follow-up analyses, Claude critically evaluates the evidence -- typically
most candidates are eliminated at this stage. A survey may end with a single
candidate worth testing, or with none."

The design commitment here is that elimination must be *cheap to justify and
expensive to avoid*. Every candidate is put through the rubric's gates, each of
which names the reason it kills. A reviewer that wants to keep a candidate has
to point at a passing gate; it cannot trade a failed gate off against an
appealing story.

Two properties this stage is built to have:

* **It can return nothing.** A stage that always produces a winner is a ranking
  function wearing a filter's clothes. If every candidate fails a gate, the
  campaign ends with none, and that is a successful run.
* **It reviews the report as well as the locus.** A report can be right about
  the locus and still assert more than the measurements carry. Overclaims are
  recorded separately from gate failures, because they say something about the
  reporting instructions rather than about the biology -- which is exactly what
  the taste stage needs to read.
"""

from __future__ import annotations

import time
from collections import Counter
from pathlib import Path

from ..agents import Harness, Session
from ..campaign import CampaignConfig, Run
from ..render import report_markdown
from ..rubric import Rubric
from ..schemas import Candidate, CandidateReport, StageRecord, Verdict

_SCHEMA_HINT = """{
  "candidate_id": "string",
  "disposition": "promote|hold|eliminate",
  "priority": 0.0,
  "failed_gates": ["gate_id"],
  "kill_reason": "one sentence, empty if not eliminated",
  "scores": {"score_id": 0.0},
  "critique": ["specific problems with the reasoning or the evidence"],
  "overclaims": ["claims asserted more strongly than the measurement supports"]
}"""

_INSTRUCTIONS = """You are reviewing a candidate report from a protein-family survey.

Your job is to eliminate. Most candidates in a survey like this one are
artefacts, rediscoveries of described systems, or associations that will not
reproduce, and the value of this review is the ones it removes. A survey that
ends with no candidate is a normal and acceptable outcome; promoting a weak
candidate is not.

How to review:
- Run every disqualifying check in the rubric below against the measurements in
  the bundle, not against the report's prose. If a check fails, eliminate and
  say which check and why.
- Do not let strong evidence on one axis excuse a failed check on another.
- Separately, compare each claim in the report against the measurement behind
  it, and record any claim asserted more strongly than that measurement
  supports. This is not grounds for elimination by itself; record it.
- Say what single observation would most change your judgement.

{rubric}
"""


def _bundle(candidate: Candidate, report: CandidateReport) -> dict:
    return {
        "candidate_id": candidate.candidate_id,
        "architecture": candidate.core_architecture,
        "lineage": candidate.lineage,
        "metrics": candidate.metrics,
        "report": {
            "title": report.title,
            "proposed_function": report.proposed_function,
            "summary": report.summary,
            "alternatives": report.alternatives,
            "n_predictions": len(report.predictions),
            "confidence": report.confidence,
        },
        "claims": [
            {"evidence_key": e.key, "statement": e.statement,
             "strength": e.strength.value, "supports": e.supports}
            for e in report.evidence
        ],
    }


def _coerce(payload: dict, candidate: Candidate, rubric: Rubric, reviewer: str) -> Verdict:
    disposition = payload.get("disposition", "eliminate")
    if disposition not in ("promote", "hold", "eliminate"):
        disposition = "eliminate"
    return Verdict(
        candidate_id=payload.get("candidate_id") or candidate.candidate_id,
        disposition=disposition,
        priority=float(payload.get("priority", 0.0) or 0.0),
        failed_gates=list(payload.get("failed_gates", [])),
        kill_reason=payload.get("kill_reason", ""),
        scores={k: float(v) for k, v in (payload.get("scores") or {}).items()},
        critique=list(payload.get("critique", [])),
        overclaims=list(payload.get("overclaims", [])),
        reviewer=reviewer,
        rubric_version=rubric.version,
    )


def run_triage(
    run: Run,
    config: CampaignConfig,
    harness: Harness,
    rubric: Rubric,
    candidates: list[Candidate],
    reports: list[CandidateReport],
) -> list[Verdict]:
    """Review every report in parallel and write the verdicts."""
    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    by_id = {c.candidate_id: c for c in candidates}
    instructions = _INSTRUCTIONS.format(rubric=rubric.as_prompt())

    pairs = [(by_id[r.candidate_id], r) for r in reports if r.candidate_id in by_id]
    sessions = [
        Session(
            session_id=f"triage.{c.candidate_id}",
            role="reviewer",
            instructions=instructions,
            task=(
                "Review this candidate. Apply every disqualifying check, decide "
                "a disposition, and record any claim the report asserts more "
                "strongly than its measurement supports."
            ),
            bundle=_bundle(c, r),
            schema_hint=_SCHEMA_HINT,
            metadata={"candidate_id": c.candidate_id},
        )
        for c, r in pairs
    ]

    results = harness.run(sessions)
    verdicts: list[Verdict] = []
    failed: list[str] = []
    for session, res in zip(sessions, results):
        cid = session.metadata["candidate_id"]
        if not res.ok:
            failed.append(f"{cid}: {res.error}")
            # A review that did not run is not an acquittal. The candidate is
            # held, never promoted, so a backend failure cannot quietly admit
            # something to the bench queue.
            verdicts.append(Verdict(
                candidate_id=cid, disposition="hold",
                kill_reason="", reviewer=f"{res.backend} (failed)",
                critique=[f"review did not complete: {res.error}"],
                rubric_version=rubric.version,
            ))
            continue
        verdicts.append(_coerce(res.output, by_id[cid], rubric, res.backend))

    # Render each report with its verdict attached.
    reports_dir = Path(run.dir) / "reports"
    verdict_by_id = {v.candidate_id: v for v in verdicts}
    banner = run.manifest.source.get("banner", "")
    for r in reports:
        md = report_markdown(r, verdict=verdict_by_id.get(r.candidate_id), banner=banner)
        (reports_dir / f"{r.candidate_id}.md").write_text(md)

    run.write_json("verdicts.json", [v.model_dump() for v in verdicts])

    counts = Counter(v.disposition for v in verdicts)
    gate_counts = Counter(g for v in verdicts for g in v.failed_gates)
    notes = [
        harness.stats.summary(),
        f"promote {counts['promote']}, hold {counts['hold']}, "
        f"eliminate {counts['eliminate']} of {len(verdicts)}",
        "most frequent gate failures: "
        + (", ".join(f"{g} x{n}" for g, n in gate_counts.most_common(5)) or "none"),
    ]
    if counts["promote"] == 0:
        notes.append(
            "no candidate survived: this is a valid outcome, not a pipeline failure"
        )
    if failed:
        notes.append(f"{len(failed)} reviews failed to run and were held, not promoted")

    run.record_stage(StageRecord(
        stage="triage",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"reports": len(reports), "rubric_version": rubric.version},
        outputs={
            "promote": counts["promote"], "hold": counts["hold"],
            "eliminate": counts["eliminate"],
            "gate_failures": dict(gate_counts),
        },
        sessions=len(sessions),
        notes=notes,
    ))
    return verdicts
