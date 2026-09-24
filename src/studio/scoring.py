"""Scoring a run against the corpus's withheld ground truth.

This is not part of the workflow -- a real campaign has no answer key. It
exists so the pipeline itself can be evaluated: does it recover the planted
novel systems, and does it kill the decoys?

The second number is the one that matters. Recovering a planted system while
also promoting the catalytically dead RT, the fragment and the diverged Cas1
would demonstrate nothing, because the article's workflow is defined by what it
discards. A run that promotes the novel systems and nothing else is the target;
a run that promotes nothing is a miss, but an honest one; a run that promotes
decoys is a failure whatever else it found.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class DecoyOutcome:
    system_id: str
    kill_reason_expected: str
    n_planted: int = 0
    n_reached_triage: int = 0
    n_killed: int = 0
    n_promoted: int = 0
    gates_that_fired: dict[str, int] = field(default_factory=dict)

    @property
    def caught(self) -> bool:
        return self.n_promoted == 0

    def line(self) -> str:
        flag = "CAUGHT" if self.caught else "LEAKED"
        gates = ", ".join(f"{g} x{n}" for g, n in
                          sorted(self.gates_that_fired.items(), key=lambda kv: -kv[1])[:3])
        return (
            f"[{flag}] {self.system_id}: {self.n_planted} planted, "
            f"{self.n_reached_triage} reached triage, {self.n_killed} killed, "
            f"{self.n_promoted} promoted"
            + (f" (by {gates})" if gates else "")
        )


@dataclass
class ScoreReport:
    run_id: str
    novel_planted: int = 0
    novel_recovered: int = 0
    novel_promoted: int = 0
    decoys: list[DecoyOutcome] = field(default_factory=list)
    described_promoted: int = 0
    unexplained_promotions: list[str] = field(default_factory=list)
    promoted_total: int = 0
    per_system_recovery: dict[str, str] = field(default_factory=dict)

    @property
    def decoys_leaked(self) -> int:
        return sum(1 for d in self.decoys if not d.caught)

    @property
    def clean(self) -> bool:
        """Every decoy killed, no described system promoted."""
        return self.decoys_leaked == 0 and self.described_promoted == 0

    def text(self) -> str:
        lines = [
            f"Scoring run {self.run_id} against withheld ground truth",
            "=" * 56,
            "",
            f"Novel systems planted   : {self.novel_planted}",
            f"  reached the candidate list: {self.novel_recovered}",
            f"  promoted to the bench     : {self.novel_promoted}",
            "",
            "Per planted system:",
        ]
        for sid, text in sorted(self.per_system_recovery.items()):
            lines.append(f"  {sid}: {text}")
        lines += ["", "Decoys:"]
        lines += ["  " + d.line() for d in self.decoys]
        lines += [
            "",
            f"Described systems promoted: {self.described_promoted} "
            f"(any number above zero is a rediscovery leak)",
            f"Total promoted            : {self.promoted_total}",
        ]
        if self.unexplained_promotions:
            lines += ["", "Promoted but not planted as novel:"]
            lines += [f"  {c}" for c in self.unexplained_promotions]
        lines += [
            "",
            f"VERDICT: {'clean' if self.clean else 'LEAKY'} "
            f"-- {self.decoys_leaked} decoy classes leaked, "
            f"{self.described_promoted} described systems promoted",
        ]
        return "\n".join(lines)


def score_run(run_dir: Path, corpus_dir: Path) -> ScoreReport:
    """Compare a run's outcome to the corpus's ``TRUTH.withheld.json``."""
    run_dir, corpus_dir = Path(run_dir), Path(corpus_dir)
    truth = json.loads((corpus_dir / "TRUTH.withheld.json").read_text())
    truth_by_anchor = {t["anchor_gene"]: t for t in truth}

    manifest = json.loads((run_dir / "manifest.json").read_text())
    candidates = json.loads((run_dir / "candidates.json").read_text()) \
        if (run_dir / "candidates.json").exists() else []
    verdicts = json.loads((run_dir / "verdicts.json").read_text()) \
        if (run_dir / "verdicts.json").exists() else []

    verdict_by_id = {v["candidate_id"]: v for v in verdicts}
    cand_by_id = {c["candidate_id"]: c for c in candidates}

    report = ScoreReport(run_id=manifest.get("run_id", run_dir.name))

    # Which planted loci made it into the candidate list, and what happened.
    reached: dict[str, list[str]] = {}
    for cid, cand in cand_by_id.items():
        anchor = cand.get("anchor_protein", "")
        t = truth_by_anchor.get(anchor)
        if t is None:
            continue
        reached.setdefault(t["system_id"], []).append(cid)

    promoted_ids = {
        cid for cid, v in verdict_by_id.items() if v.get("disposition") == "promote"
    }
    report.promoted_total = len(promoted_ids)

    # Novel systems.
    novel_truth = [t for t in truth if t["is_novel"]]
    report.novel_planted = len(novel_truth)
    novel_ids = {t["system_id"] for t in novel_truth}
    for sid in sorted(novel_ids):
        planted = sum(1 for t in novel_truth if t["system_id"] == sid)
        got = reached.get(sid, [])
        prom = [c for c in got if c in promoted_ids]
        report.novel_recovered += len(got)
        report.novel_promoted += len(prom)
        report.per_system_recovery[sid] = (
            f"{planted} planted, {len(got)} became candidates, {len(prom)} promoted"
        )

    # Decoys.
    decoy_ids = sorted({t["system_id"] for t in truth if t["system_id"].startswith("decoy:")})
    for sid in decoy_ids:
        planted = [t for t in truth if t["system_id"] == sid]
        got = reached.get(sid, [])
        prom = [c for c in got if c in promoted_ids]
        gates: Counter = Counter()
        for cid in got:
            for g in verdict_by_id.get(cid, {}).get("failed_gates", []):
                gates[g] += 1
        report.decoys.append(DecoyOutcome(
            system_id=sid,
            kill_reason_expected=planted[0].get("kill_reason", ""),
            n_planted=len(planted),
            n_reached_triage=len(got),
            n_killed=len(got) - len(prom),
            n_promoted=len(prom),
            gates_that_fired=dict(gates),
        ))

    # Described systems promoted: a rediscovery leak.
    for cid in promoted_ids:
        anchor = cand_by_id.get(cid, {}).get("anchor_protein", "")
        t = truth_by_anchor.get(anchor)
        if t is None:
            report.unexplained_promotions.append(f"{cid} (anchor not a planted locus)")
        elif not t["is_novel"] and not t["system_id"].startswith("decoy:"):
            report.described_promoted += 1
            report.unexplained_promotions.append(f"{cid} ({t['system_id']})")

    return report
