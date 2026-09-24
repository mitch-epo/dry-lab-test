"""Stage 7 -- taste: study the hypotheses, and rewrite the instructions.

"Because Claude produces hypotheses so prolifically, the hypotheses themselves
have become an object of study for us. With hundreds to thousands of candidate
reports from a single campaign, we have been asking what distinguishes the
proposals we judge worth testing from those we set aside. What we learn goes
back into the instructions we give Claude."

This stage mines the campaign's own report corpus for features that separate
promoted from set-aside candidates, and emits the next rubric version.

THE CIRCULARITY PROBLEM, AND WHAT IS DONE ABOUT IT
--------------------------------------------------
The verdicts were produced *by* the rubric. Learning "what distinguishes
promoted from set-aside" from rubric-generated verdicts and then feeding that
back into the rubric is a loop that re-derives its own premises: features the
rubric already weights will always look discriminative, and tightening them
will always look justified. Run unexamined, this manufactures confidence.

So the stage separates its findings by whether they can survive that critique:

* **External labels take precedence.** Outcomes recorded by people -- bench
  results, or a reviewer overriding a verdict -- are read from ``labels.yaml``
  when present. Analysis against external labels is not circular and is
  reported as ``external``.

* **Without external labels, the stage restricts what it will act on.** It
  proposes changes only where the rubric-generated verdicts still carry
  information about the rubric itself rather than about the biology:
  - features that discriminate but are *not* in the rubric (the rubric is
    responding to something it does not name, which is worth naming);
  - gates that never fire, or fire on everything (a gate that eliminates
    nothing is not doing work; one that eliminates everything is mis-set);
  - systematic overclaims in the reports, which are a property of the
    reporting instructions and are measured independently of the verdicts.
  Weight changes to criteria already in the rubric are computed and shown, but
  marked ``circular`` and not applied.

Every proposal records which category it came from, and the emitted rubric
records why each change was made.
"""

from __future__ import annotations

import json
import math
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from ..campaign import CampaignConfig, Run
from ..rubric import Rubric, Score
from ..schemas import CandidateReport, StageRecord, Verdict


@dataclass
class FeatureDiscrimination:
    """How well one measured feature separates promoted from set-aside."""

    feature: str
    auroc: float                 # 0.5 = no signal, 1.0 = perfect separation
    n_promoted: int
    n_set_aside: int
    mean_promoted: float
    mean_set_aside: float
    in_rubric: bool
    kind: str = "numeric"

    @property
    def strength(self) -> float:
        """Distance from chance, in 0..0.5."""
        return abs(self.auroc - 0.5)

    def line(self) -> str:
        flag = "in rubric" if self.in_rubric else "NOT in rubric"
        return (
            f"{self.feature}: AUROC {self.auroc:.2f} ({flag}); "
            f"promoted mean {self.mean_promoted:.3g} vs "
            f"set-aside {self.mean_set_aside:.3g}"
        )


@dataclass
class Proposal:
    """A change to the rubric, with its justification and whether it was applied."""

    kind: str                    # add_score | reweight | retire_gate | loosen_gate | guidance
    target: str
    detail: str
    evidence: str
    basis: str                   # external | structural | circular
    applied: bool = False


@dataclass
class TasteResult:
    n_reports: int = 0
    n_promoted: int = 0
    n_set_aside: int = 0
    label_source: str = "verdicts"
    discriminations: list[dict] = field(default_factory=list)
    gate_activity: dict[str, int] = field(default_factory=dict)
    inert_gates: list[str] = field(default_factory=list)
    saturated_gates: list[str] = field(default_factory=list)
    overclaim_patterns: dict[str, int] = field(default_factory=dict)
    proposals: list[dict] = field(default_factory=list)
    next_rubric_version: int = 0
    notes: list[str] = field(default_factory=list)

    def report(self) -> str:
        lines = [
            "Taste analysis",
            "=" * 14,
            "",
            f"corpus: {self.n_reports} reports, {self.n_promoted} promoted, "
            f"{self.n_set_aside} set aside",
            f"labels: {self.label_source}",
            "",
            "Feature discrimination (ranked by distance from chance):",
        ]
        for d in self.discriminations[:12]:
            lines.append("  " + FeatureDiscrimination(**d).line())
        lines += ["", "Gate activity:"]
        for g, n in sorted(self.gate_activity.items(), key=lambda kv: -kv[1]):
            lines.append(f"  {g}: eliminated {n}")
        if self.inert_gates:
            lines.append(f"  inert (never fired): {', '.join(self.inert_gates)}")
        if self.saturated_gates:
            lines.append(f"  saturated (fired on everything): {', '.join(self.saturated_gates)}")
        if self.overclaim_patterns:
            lines += ["", "Recurring overclaims:"]
            for k, n in sorted(self.overclaim_patterns.items(), key=lambda kv: -kv[1]):
                lines.append(f"  {k}: {n}")
        lines += ["", "Proposals:"]
        for p in self.proposals:
            mark = "APPLIED" if p["applied"] else "held"
            lines.append(f"  [{mark}/{p['basis']}] {p['kind']} {p['target']}: {p['detail']}")
        return "\n".join(lines)


def auroc(positive: list[float], negative: list[float]) -> float:
    """Area under the ROC curve via the Mann-Whitney U statistic.

    Ties count a half, which matters here: many of these features are boolean
    or heavily quantised, and ignoring ties would inflate the separation.
    """
    if not positive or not negative:
        return 0.5
    wins = 0.0
    for p in positive:
        for n in negative:
            if p > n:
                wins += 1.0
            elif p == n:
                wins += 0.5
    return wins / (len(positive) * len(negative))


def _numeric(value: Any) -> float | None:
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _rubric_features(rubric: Rubric) -> set[str]:
    """Metric names the rubric already refers to, read out of its expressions."""
    import ast as _ast

    names: set[str] = set()
    for expr in [g.expr for g in rubric.gates] + [s.expr for s in rubric.scores]:
        try:
            tree = _ast.parse(expr, mode="eval")
        except SyntaxError:
            continue
        for node in _ast.walk(tree):
            if isinstance(node, _ast.Name):
                names.add(node.id)
    return names


def load_labels(path: Path) -> dict[str, str]:
    """External outcome labels: ``{candidate_id: promote|set_aside}``.

    Written by people -- a bench result, or a reviewer overriding a verdict.
    When this file exists the analysis is not circular, and the stage says so.
    """
    if not Path(path).exists():
        return {}
    doc = yaml.safe_load(Path(path).read_text()) or {}
    out: dict[str, str] = {}
    for cid, value in (doc.get("labels") or {}).items():
        v = str(value).lower()
        if v in ("promote", "promoted", "tested", "worth_testing", "true", "yes"):
            out[cid] = "promote"
        elif v in ("set_aside", "eliminate", "eliminated", "false", "no"):
            out[cid] = "set_aside"
    return out


def analyse_corpus(bundle: dict[str, Any], rubric: Rubric) -> dict[str, Any]:
    """The analyst role: rank features by how well they separate the two classes.

    Called through the harness so the meta-analysis is itself a recorded
    session, with the same provenance as everything else in the campaign.
    """
    rows = bundle.get("rows", [])
    known = _rubric_features(rubric)

    by_feature: dict[str, tuple[list[float], list[float]]] = defaultdict(lambda: ([], []))
    n_pos = n_neg = 0
    for row in rows:
        promoted = bool(row.get("promoted"))
        n_pos += promoted
        n_neg += not promoted
        for key, raw in (row.get("metrics") or {}).items():
            v = _numeric(raw)
            if v is None:
                continue
            by_feature[key][0 if promoted else 1].append(v)

    discriminations: list[FeatureDiscrimination] = []
    for feature, (pos, neg) in by_feature.items():
        if len(pos) < 2 or len(neg) < 2:
            continue
        a = auroc(pos, neg)
        discriminations.append(FeatureDiscrimination(
            feature=feature,
            auroc=round(a, 4),
            n_promoted=len(pos),
            n_set_aside=len(neg),
            mean_promoted=round(sum(pos) / len(pos), 4),
            mean_set_aside=round(sum(neg) / len(neg), 4),
            in_rubric=feature in known,
            kind="boolean" if set(pos + neg) <= {0.0, 1.0} else "numeric",
        ))
    discriminations.sort(key=lambda d: -d.strength)
    return {
        "n_promoted": n_pos,
        "n_set_aside": n_neg,
        "discriminations": [asdict(d) for d in discriminations],
    }


def _collect_rows(run: Run, labels: dict[str, str]) -> tuple[list[dict], str]:
    """Build the analysis table from this run's candidates and verdicts."""
    candidates = {c["candidate_id"]: c for c in (run.read_json("candidates.json") or [])}
    verdicts = {v["candidate_id"]: v for v in (run.read_json("verdicts.json") or [])}

    rows: list[dict] = []
    used_external = False
    for cid, cand in candidates.items():
        v = verdicts.get(cid)
        if v is None:
            continue
        if cid in labels:
            promoted = labels[cid] == "promote"
            used_external = True
        else:
            promoted = v["disposition"] == "promote"
        rows.append({
            "candidate_id": cid,
            "promoted": promoted,
            "metrics": cand.get("metrics", {}),
            "failed_gates": v.get("failed_gates", []),
            "overclaims": v.get("overclaims", []),
        })
    return rows, ("external" if used_external else "verdicts")


def run_taste(
    run: Run,
    config: CampaignConfig,
    harness,
    rubric: Rubric,
    *,
    labels_path: Path | None = None,
    min_auroc: float = 0.70,
    min_corpus: int = 8,
) -> TasteResult:
    """Mine the report corpus and emit the next rubric version."""
    from ..agents import Session

    t0 = time.monotonic()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    labels = load_labels(labels_path) if labels_path else {}
    rows, label_source = _collect_rows(run, labels)

    result = TasteResult(
        n_reports=len(rows),
        n_promoted=sum(1 for r in rows if r["promoted"]),
        n_set_aside=sum(1 for r in rows if not r["promoted"]),
        label_source=label_source,
    )

    session = Session(
        session_id="taste.corpus-analysis",
        role="analyst",
        instructions=(
            "You are analysing a campaign's own report corpus to find what "
            "separates the proposals judged worth testing from those set "
            "aside. Report feature separation, not a story about it. Note "
            "explicitly where the labels were produced by the same rubric you "
            "are being asked to revise.\n\n" + rubric.as_prompt()
        ),
        task=(
            "Rank the measured features by how well they separate promoted "
            "from set-aside candidates."
        ),
        bundle={"rows": rows},
        max_tokens=8192,
    )
    analysis = harness.run_one(session)
    payload = analysis.output if analysis.ok else analyse_corpus(session.bundle, rubric)
    result.discriminations = payload.get("discriminations", [])

    # -- gate activity -------------------------------------------------------
    gate_counts = Counter(g for r in rows for g in r["failed_gates"])
    result.gate_activity = dict(gate_counts)
    all_gate_ids = [g.gate_id for g in rubric.gates if g.enabled]
    result.inert_gates = [g for g in all_gate_ids if gate_counts.get(g, 0) == 0]
    result.saturated_gates = [
        g for g in all_gate_ids if rows and gate_counts.get(g, 0) == len(rows)
    ]

    # -- overclaim patterns --------------------------------------------------
    patterns: Counter = Counter()
    for r in rows:
        for text in r["overclaims"]:
            # Overclaims read "<key>: asserted as X but ...". Group on the key.
            key = text.split(":", 1)[0].strip() if ":" in text else text
            patterns[key[:80]] += 1
    result.overclaim_patterns = dict(patterns)

    # -- proposals -----------------------------------------------------------
    proposals: list[Proposal] = []
    next_rubric = Rubric.from_dict(rubric.to_dict())
    external = label_source == "external"
    enough = len(rows) >= min_corpus

    for d in result.discriminations:
        fd = FeatureDiscrimination(**d)
        if fd.auroc < min_auroc and fd.auroc > 1 - min_auroc:
            continue
        if fd.in_rubric:
            proposals.append(Proposal(
                kind="reweight",
                target=fd.feature,
                detail=(
                    f"separation AUROC {fd.auroc:.2f} suggests raising this "
                    f"criterion's weight"
                ),
                evidence=fd.line(),
                basis="external" if external else "circular",
                applied=False,
            ))
            continue
        # A feature outside the rubric that separates the classes means the
        # reviewer is responding to something the rubric does not name. Naming
        # it is a real improvement even when the labels came from the rubric.
        direction = "" if fd.auroc >= 0.5 else " (inverted: lower is better)"
        expr = (
            f"clamp({fd.feature})" if 0.0 <= fd.mean_promoted <= 1.0 and fd.kind != "boolean"
            else (fd.feature if fd.kind == "boolean" else f"clamp({fd.feature} / "
                  f"{max(abs(fd.mean_promoted), 1e-6):.6g})")
        )
        if fd.auroc < 0.5:
            expr = f"1.0 - ({expr})"
        proposals.append(Proposal(
            kind="add_score",
            target=fd.feature,
            detail=f"add as a ranking criterion with weight 0.5{direction}",
            evidence=fd.line(),
            basis="external" if external else "structural",
            applied=enough,
        ))
        if enough:
            next_rubric.scores.append(Score(
                score_id=f"learned_{fd.feature}",
                question=(
                    f"Does {fd.feature} favour this candidate? "
                    f"(learned from the campaign corpus)"
                ),
                expr=expr,
                weight=0.5,
                rationale=(
                    f"Separated promoted from set-aside at AUROC {fd.auroc:.2f} "
                    f"over {fd.n_promoted + fd.n_set_aside} reports in rubric "
                    f"v{rubric.version}, and was not named by that rubric."
                ),
            ))

    for gate_id in result.inert_gates:
        proposals.append(Proposal(
            kind="retire_gate",
            target=gate_id,
            detail=(
                "eliminated nothing in this corpus; either it is redundant with "
                "another gate or its threshold is unreachable"
            ),
            evidence=f"0 of {len(rows)} candidates failed it",
            basis="structural",
            applied=False,   # never retired automatically: a quiet gate may be
                             # quiet because upstream filtering already works.
        ))
    for gate_id in result.saturated_gates:
        proposals.append(Proposal(
            kind="loosen_gate",
            target=gate_id,
            detail="eliminated every candidate; the threshold cannot be met as written",
            evidence=f"{len(rows)} of {len(rows)} candidates failed it",
            basis="structural",
            applied=False,
        ))

    for pattern, n in patterns.most_common(3):
        if n < max(2, len(rows) // 10):
            continue
        guidance = (
            f"Rate the {pattern!r} evidence from its measurement, not from the "
            f"strength of the overall case: it was over-stated in {n} reports "
            f"under rubric v{rubric.version}."
        )
        proposals.append(Proposal(
            kind="guidance",
            target=pattern,
            detail=guidance,
            evidence=f"{n} reports overclaimed this",
            basis="structural",
            applied=True,
        ))
        if guidance not in next_rubric.guidance:
            next_rubric.guidance.append(guidance)

    applied = [p for p in proposals if p.applied]
    if applied:
        next_rubric.version = rubric.version + 1
        next_rubric.history = list(rubric.history) + [{
            "from_version": rubric.version,
            "run_id": run.run_id,
            "label_source": label_source,
            "corpus_size": len(rows),
            "applied": [
                {"kind": p.kind, "target": p.target, "basis": p.basis,
                 "evidence": p.evidence}
                for p in applied
            ],
        }]
        next_rubric.notes = list(rubric.notes) + [
            f"v{next_rubric.version} derived from run {run.run_id} "
            f"({len(rows)} reports, labels: {label_source})."
        ]
        next_rubric.save(Path(run.dir) / "rubric.next.yaml")
        result.next_rubric_version = next_rubric.version
    else:
        result.next_rubric_version = rubric.version

    result.proposals = [asdict(p) for p in proposals]
    result.notes = [
        f"{len(rows)} reports analysed; labels from {label_source}",
        f"{len(applied)} of {len(proposals)} proposals applied",
    ]
    if not external:
        result.notes.append(
            "labels came from the same rubric being revised, so weight changes "
            "to existing criteria are reported but NOT applied; supply "
            "labels.yaml with external outcomes to lift that restriction"
        )
    if not enough:
        result.notes.append(
            f"corpus of {len(rows)} is below the minimum of {min_corpus}: "
            f"feature proposals are reported but not applied"
        )
    if result.n_promoted == 0:
        result.notes.append(
            "nothing was promoted, so no contrast exists to learn from; "
            "discrimination figures are undefined"
        )

    run.write_json("taste.json", asdict(result))
    run.write_text("taste.txt", result.report() + "\n")
    run.record_stage(StageRecord(
        stage="taste",
        started_at=started,
        finished_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        duration_s=time.monotonic() - t0,
        inputs={"reports": len(rows), "rubric_version": rubric.version,
                "label_source": label_source},
        outputs={"proposals": len(proposals), "applied": len(applied),
                 "next_rubric_version": result.next_rubric_version},
        sessions=1,
        notes=result.notes,
    ))
    return result
