"""The campaign data model.

The article's workflow turns on one discipline: a candidate report "proposes a
function and describes the evidence supporting its claims", and a later pass
"critically evaluates the evidence". That only works if evidence is a first-
class object rather than prose. So every claim here carries the measurement it
rests on, the method that produced it, and a handle to the record it came
from -- which is what makes the triage stage able to disagree with the report
stage on the facts rather than on the writing.
"""

from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class Strength(str, Enum):
    """How much weight a single piece of evidence can carry."""

    DIRECT = "direct"            # a measurement of the thing claimed
    STRONG = "strong"            # a well-powered statistical association
    SUGGESTIVE = "suggestive"    # consistent with, but also with alternatives
    WEAK = "weak"                # anecdotal or underpowered
    ABSENT = "absent"            # the check was run and found nothing


class Evidence(BaseModel):
    """One measurement supporting or undercutting a claim."""

    model_config = ConfigDict(extra="forbid")

    key: str                                  # stable identifier, e.g. "catalytic_core"
    statement: str                            # human-readable, quotable in a report
    value: Any = None                         # the number or flag behind the statement
    strength: Strength = Strength.SUGGESTIVE
    method: str = ""                          # how it was computed
    records: list[str] = Field(default_factory=list)   # accessions / gene ids
    supports: bool = True                     # False when it argues against the claim

    def line(self) -> str:
        mark = "+" if self.supports else "-"
        rec = f" [{', '.join(self.records[:4])}]" if self.records else ""
        return f"{mark} ({self.strength.value}) {self.statement}{rec}"


class Candidate(BaseModel):
    """A locus that the search stage judged worth writing up."""

    model_config = ConfigDict(extra="forbid")

    candidate_id: str
    anchor_protein: str
    contig_id: str
    lineage: str = ""
    architecture: str = ""
    core_architecture: str = ""
    group_signature: str = ""
    novelty: float = 0.0
    closest_system: str = ""
    closest_system_score: float = 0.0
    n_family_members: int = 0
    metrics: dict[str, Any] = Field(default_factory=dict)
    evidence: list[Evidence] = Field(default_factory=list)
    homolog_loci: list[str] = Field(default_factory=list)

    def metric(self, key: str, default: Any = None) -> Any:
        return self.metrics.get(key, default)


class Prediction(BaseModel):
    """A falsifiable consequence of the proposed function."""

    model_config = ConfigDict(extra="forbid")

    statement: str
    assay: str
    would_falsify: str


class CandidateReport(BaseModel):
    """The short, human-readable report the workflow writes per candidate."""

    model_config = ConfigDict(extra="forbid")

    candidate_id: str
    title: str
    proposed_function: str
    summary: str
    evidence: list[Evidence] = Field(default_factory=list)
    alternatives: list[str] = Field(default_factory=list)
    predictions: list[Prediction] = Field(default_factory=list)
    open_questions: list[str] = Field(default_factory=list)
    confidence: float = 0.0
    author: str = ""
    rubric_version: int = 0
    created_at: str = Field(default_factory=_now)

    @property
    def supporting(self) -> list[Evidence]:
        return [e for e in self.evidence if e.supports]

    @property
    def contradicting(self) -> list[Evidence]:
        return [e for e in self.evidence if not e.supports]


Disposition = Literal["promote", "hold", "eliminate"]


class Verdict(BaseModel):
    """The triage stage's judgement on one report."""

    model_config = ConfigDict(extra="forbid")

    candidate_id: str
    disposition: Disposition
    priority: float = 0.0
    failed_gates: list[str] = Field(default_factory=list)
    kill_reason: str = ""
    scores: dict[str, float] = Field(default_factory=dict)
    critique: list[str] = Field(default_factory=list)
    overclaims: list[str] = Field(default_factory=list)
    reviewer: str = ""
    rubric_version: int = 0
    created_at: str = Field(default_factory=_now)

    @property
    def survived(self) -> bool:
        return self.disposition == "promote"


class ReproductionResult(BaseModel):
    """One positive control from the methods check."""

    model_config = ConfigDict(extra="forbid")

    system_id: str
    expected: int
    recovered: int
    false_positives: int
    recall: float
    precision: float
    min_recall: float
    min_precision: float

    @property
    def passed(self) -> bool:
        return self.recall >= self.min_recall and self.precision >= self.min_precision

    def line(self) -> str:
        flag = "PASS" if self.passed else "FAIL"
        return (
            f"[{flag}] {self.system_id}: recall {self.recall:.2f} "
            f"(>= {self.min_recall:.2f}), precision {self.precision:.2f} "
            f"(>= {self.min_precision:.2f}) on {self.expected} instances"
        )


class StageRecord(BaseModel):
    """One executed stage, for the run manifest."""

    model_config = ConfigDict(extra="forbid")

    stage: str
    started_at: str
    finished_at: str = ""
    duration_s: float = 0.0
    inputs: dict[str, Any] = Field(default_factory=dict)
    outputs: dict[str, Any] = Field(default_factory=dict)
    sessions: int = 0
    notes: list[str] = Field(default_factory=list)


class Manifest(BaseModel):
    """Everything needed to say what a run did and what it rests on."""

    model_config = ConfigDict(extra="forbid")

    campaign: str
    run_id: str
    created_at: str = Field(default_factory=_now)
    source: dict[str, Any] = Field(default_factory=dict)
    synthetic_data: bool = False
    rubric_version: int = 0
    backend: str = ""
    stages: list[StageRecord] = Field(default_factory=list)
    gate_passed: bool | None = None
    counts: dict[str, int] = Field(default_factory=dict)

    def stage(self, name: str) -> StageRecord | None:
        return next((s for s in self.stages if s.stage == name), None)


class BenchPackage(BaseModel):
    """The wet-lab handoff for a candidate that survived triage."""

    model_config = ConfigDict(extra="forbid")

    candidate_id: str
    title: str
    proposed_function: str
    architecture: str = ""
    supporting_loci: list[str] = Field(default_factory=list)
    n_instances: int = 1
    constructs: list[dict[str, str]] = Field(default_factory=list)
    assays: list[dict[str, str]] = Field(default_factory=list)
    controls: list[str] = Field(default_factory=list)
    falsification: list[str] = Field(default_factory=list)
    biosafety_notes: list[str] = Field(default_factory=list)
    created_at: str = Field(default_factory=_now)
