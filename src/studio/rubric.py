"""The rubric: the studio's scientific taste, written down and versioned.

This is the pivot of the whole design. The article ends on it: hundreds to
thousands of candidate reports per campaign, a study of what separates the
proposals worth testing from the ones set aside, and what is learned going
"back into the instructions we give Claude". For that loop to close, the
instructions have to be an artifact the pipeline can both *read* and *rewrite*
-- not prose buried in a prompt string.

So a rubric is data. It has two parts:

* **gates** -- boolean conditions a candidate must satisfy. A failed gate
  eliminates the candidate and names the reason. This is where "typically most
  candidates are eliminated at this stage" actually happens.
* **scores** -- expressions in 0..1 with weights, combined into the priority
  that ranks survivors.

Both are expressions over a candidate's measured evidence, evaluated through a
restricted AST walker -- no ``eval``. The same rubric renders to the prompt
text an LLM reviewer receives, so the deterministic and the model-backed
reviewers are held to one standard, and :mod:`studio.stages.taste` can revise
that standard from the campaign's own history.
"""

from __future__ import annotations

import ast
import operator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# --- a restricted expression evaluator ---------------------------------------

_BIN_OPS = {
    ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
    ast.Div: operator.truediv, ast.Pow: operator.pow, ast.Mod: operator.mod,
}
_CMP_OPS = {
    ast.Eq: operator.eq, ast.NotEq: operator.ne, ast.Lt: operator.lt,
    ast.LtE: operator.le, ast.Gt: operator.gt, ast.GtE: operator.ge,
    ast.In: lambda a, b: a in b, ast.NotIn: lambda a, b: a not in b,
}
_FUNCS: dict[str, Any] = {
    "min": min, "max": max, "abs": abs, "len": len, "round": round,
    "float": float, "int": int, "bool": bool, "any": any, "all": all,
    "clamp": lambda x, lo=0.0, hi=1.0: max(lo, min(hi, x)),
}


class RubricExpressionError(ValueError):
    """Raised when an expression is malformed or uses an unavailable name."""


def evaluate(expr: str, variables: dict[str, Any]) -> Any:
    """Evaluate a rubric expression against a variable bag.

    Only literals, names from ``variables``, the arithmetic and comparison
    operators above, boolean logic, conditionals and the whitelisted functions
    are permitted. Anything else -- attribute access, subscripting a callable,
    imports, comprehensions over arbitrary objects -- raises. The rubric is a
    file the taste stage rewrites automatically, so it must never be a path to
    arbitrary code execution.
    """
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as exc:
        raise RubricExpressionError(f"cannot parse {expr!r}: {exc}") from exc

    def walk(node: ast.AST) -> Any:
        if isinstance(node, ast.Expression):
            return walk(node.body)
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            if node.id in variables:
                return variables[node.id]
            if node.id in _FUNCS:
                return _FUNCS[node.id]
            raise RubricExpressionError(f"unknown name {node.id!r} in {expr!r}")
        if isinstance(node, ast.UnaryOp):
            v = walk(node.operand)
            if isinstance(node.op, ast.Not):
                return not v
            if isinstance(node.op, ast.USub):
                return -v
            if isinstance(node.op, ast.UAdd):
                return +v
            raise RubricExpressionError(f"unsupported unary op in {expr!r}")
        if isinstance(node, ast.BinOp):
            fn = _BIN_OPS.get(type(node.op))
            if fn is None:
                raise RubricExpressionError(f"unsupported operator in {expr!r}")
            return fn(walk(node.left), walk(node.right))
        if isinstance(node, ast.BoolOp):
            values = [walk(v) for v in node.values]
            return all(values) if isinstance(node.op, ast.And) else any(values)
        if isinstance(node, ast.Compare):
            left = walk(node.left)
            for op, comparator in zip(node.ops, node.comparators):
                fn = _CMP_OPS.get(type(op))
                if fn is None:
                    raise RubricExpressionError(f"unsupported comparison in {expr!r}")
                right = walk(comparator)
                if not fn(left, right):
                    return False
                left = right
            return True
        if isinstance(node, ast.IfExp):
            return walk(node.body) if walk(node.test) else walk(node.orelse)
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in _FUNCS:
                raise RubricExpressionError(f"only whitelisted functions may be called in {expr!r}")
            return _FUNCS[node.func.id](*[walk(a) for a in node.args])
        if isinstance(node, (ast.List, ast.Tuple)):
            return [walk(e) for e in node.elts]
        raise RubricExpressionError(f"unsupported syntax {type(node).__name__} in {expr!r}")

    return walk(tree)


# --- rubric objects ----------------------------------------------------------


@dataclass
class Gate:
    """A condition a candidate must satisfy to survive triage."""

    gate_id: str
    question: str
    expr: str
    kill_reason: str
    rationale: str = ""
    enabled: bool = True

    def check(self, variables: dict[str, Any]) -> tuple[bool, str]:
        # A gate that cannot be evaluated is a fault in the pipeline, not in the
        # candidate, so RubricExpressionError propagates rather than being
        # caught and reported as a kill. Swallowing it here would silently
        # discard real candidates whenever a metric went missing.
        ok = bool(evaluate(self.expr, variables))
        return ok, ("" if ok else self.kill_reason)


@dataclass
class Score:
    """A weighted 0..1 contribution to a candidate's priority."""

    score_id: str
    question: str
    expr: str
    weight: float = 1.0
    rationale: str = ""

    def value(self, variables: dict[str, Any]) -> float:
        v = evaluate(self.expr, variables)
        if isinstance(v, bool):
            v = 1.0 if v else 0.0
        return max(0.0, min(1.0, float(v)))


@dataclass
class Rubric:
    """A versioned set of gates, scores and written guidance."""

    version: int = 1
    name: str = "default"
    notes: list[str] = field(default_factory=list)
    gates: list[Gate] = field(default_factory=list)
    scores: list[Score] = field(default_factory=list)
    promote_threshold: float = 0.55
    hold_threshold: float = 0.35
    guidance: list[str] = field(default_factory=list)
    history: list[dict] = field(default_factory=list)

    # -- evaluation ---------------------------------------------------------

    def apply(self, variables: dict[str, Any]) -> tuple[list[str], list[str], dict[str, float], float]:
        """Run the rubric. Returns ``(failed_gate_ids, reasons, scores, priority)``."""
        failed: list[str] = []
        reasons: list[str] = []
        for gate in self.gates:
            if not gate.enabled:
                continue
            ok, reason = gate.check(variables)
            if not ok:
                failed.append(gate.gate_id)
                reasons.append(reason)

        scores: dict[str, float] = {}
        total_w = 0.0
        acc = 0.0
        for s in self.scores:
            v = s.value(variables)
            scores[s.score_id] = v
            acc += v * s.weight
            total_w += s.weight
        priority = acc / total_w if total_w else 0.0
        return failed, reasons, scores, priority

    def disposition(self, failed: list[str], priority: float) -> str:
        if failed:
            return "eliminate"
        if priority >= self.promote_threshold:
            return "promote"
        if priority >= self.hold_threshold:
            return "hold"
        return "eliminate"

    # -- rendering for a model-backed reviewer ------------------------------

    def as_prompt(self) -> str:
        """Render the rubric as the instruction block an LLM reviewer receives.

        The deterministic backend evaluates the expressions; a model backend
        reads this text. They are generated from one object so the two cannot
        drift apart, and a taste-stage revision reaches both at once.
        """
        lines = [
            f"# Review rubric v{self.version} ({self.name})",
            "",
            "## Disqualifying checks",
            "Run every check. If any fails, the candidate is eliminated and you "
            "state which check failed and why. Do not trade a failed check off "
            "against strong evidence elsewhere.",
            "",
        ]
        for g in self.gates:
            if not g.enabled:
                continue
            lines.append(f"- **{g.gate_id}** -- {g.question}")
            lines.append(f"  - eliminate if: `{g.expr}` is false")
            lines.append(f"  - reason to give: {g.kill_reason}")
            if g.rationale:
                lines.append(f"  - why this matters: {g.rationale}")
        lines += ["", "## Ranking criteria", ""]
        for s in self.scores:
            lines.append(f"- **{s.score_id}** (weight {s.weight:g}) -- {s.question}")
            lines.append(f"  - measured as: `{s.expr}`")
            if s.rationale:
                lines.append(f"  - why this matters: {s.rationale}")
        lines += [
            "",
            f"Promote at priority >= {self.promote_threshold:.2f}; "
            f"hold between {self.hold_threshold:.2f} and that; otherwise eliminate.",
        ]
        if self.guidance:
            lines += ["", "## Learned guidance", ""]
            lines += [f"- {g}" for g in self.guidance]
        return "\n".join(lines)

    # -- persistence --------------------------------------------------------

    def to_dict(self) -> dict:
        return {
            "version": self.version,
            "name": self.name,
            "notes": self.notes,
            "promote_threshold": self.promote_threshold,
            "hold_threshold": self.hold_threshold,
            "gates": [
                {"gate_id": g.gate_id, "question": g.question, "expr": g.expr,
                 "kill_reason": g.kill_reason, "rationale": g.rationale,
                 "enabled": g.enabled}
                for g in self.gates
            ],
            "scores": [
                {"score_id": s.score_id, "question": s.question, "expr": s.expr,
                 "weight": s.weight, "rationale": s.rationale}
                for s in self.scores
            ],
            "guidance": self.guidance,
            "history": self.history,
        }

    def save(self, path: Path) -> None:
        Path(path).write_text(yaml.safe_dump(self.to_dict(), sort_keys=False))

    @classmethod
    def from_dict(cls, d: dict) -> "Rubric":
        return cls(
            version=int(d.get("version", 1)),
            name=d.get("name", "default"),
            notes=list(d.get("notes", [])),
            promote_threshold=float(d.get("promote_threshold", 0.55)),
            hold_threshold=float(d.get("hold_threshold", 0.35)),
            gates=[Gate(**g) for g in d.get("gates", [])],
            scores=[Score(**s) for s in d.get("scores", [])],
            guidance=list(d.get("guidance", [])),
            history=list(d.get("history", [])),
        )

    @classmethod
    def load(cls, path: Path) -> "Rubric":
        return cls.from_dict(yaml.safe_load(Path(path).read_text()))


def default_rubric() -> Rubric:
    """The rubric a campaign starts from, before the taste loop has run.

    Every gate here corresponds to a way a comparative-genomics hypothesis
    routinely fails. They are stated as checks rather than as advice because
    advice is not auditable: a run has to be able to show which check killed
    which candidate.
    """
    return Rubric(
        version=1,
        name="rt-campaign-baseline",
        notes=[
            "Baseline rubric. Gates encode the failure modes listed in the "
            "literature digest; weights are uncalibrated until the taste stage "
            "has a labelled corpus of reports to learn from.",
        ],
        gates=[
            Gate(
                gate_id="catalytic_core",
                question="Is the enzyme's catalytic core intact?",
                expr="triad_intact",
                kill_reason=(
                    "the catalytic triad is disrupted, so the protein cannot "
                    "perform the chemistry the proposed function requires"
                ),
                rationale=(
                    "A profile hit says a protein belongs to a family; it does "
                    "not say the protein still works. Pseudogenes and dead "
                    "paralogues score well against the profile."
                ),
            ),
            Gate(
                gate_id="domain_coverage",
                question="Does the hit cover the domain's conserved core?",
                expr="core_coverage >= 0.7 and aligned_fraction >= 0.35",
                kill_reason=(
                    "the hit covers only part of the domain: this is a fragment "
                    "or a mis-called gene boundary, not a family member"
                ),
                rationale=(
                    "A partial hit to a long domain is the commonest false "
                    "positive. Coverage of the profile's high-information "
                    "columns is what settles it; the fraction of total profile "
                    "length is only a floor against absurdly short alignments. "
                    "Weighting the two equally punishes divergence rather than "
                    "incompleteness -- a genuinely diverged member aligns "
                    "across the conserved core and not across the variable "
                    "flanks, which is what being diverged means."
                ),
            ),
            Gate(
                gate_id="not_low_complexity",
                question="Is the hit driven by genuine homology rather than composition?",
                expr="low_complexity <= 0.3",
                kill_reason=(
                    "the protein is largely low-complexity, so the profile score "
                    "reflects amino-acid composition rather than homology"
                ),
            ),
            Gate(
                gate_id="significance",
                question="Is the family assignment statistically sound?",
                expr="evalue <= 1e-5 and profile_z >= 8",
                kill_reason="the family assignment is not significant at the campaign threshold",
            ),
            Gate(
                gate_id="novelty",
                question="Does the locus fit a described system?",
                expr="novelty >= 0.35",
                kill_reason=(
                    "the locus is accounted for by a described system, so there "
                    "is nothing new to propose"
                ),
                rationale=(
                    "Rediscovery is the default outcome of a family survey and "
                    "has to be ruled out explicitly, not assumed away."
                ),
            ),
            Gate(
                gate_id="partner_is_not_known",
                question="Is the 'unknown' partner actually a diverged known protein?",
                expr="partner_known_neglog_e < 3.0",
                kill_reason=(
                    "the partner gene is a diverged member of a described "
                    "family, so the locus is a known system with a weak "
                    "annotation, not a new one"
                ),
                rationale=(
                    "A strict label threshold silently converts diverged "
                    "members of described families into apparent unknowns. This "
                    "is the single most productive way to fool an architecture "
                    "comparison. Measured as -log10(E) so it is comparable "
                    "across partner families of different lengths."
                ),
            ),
            Gate(
                gate_id="reproducible_association",
                question="Is the association seen more than once?",
                expr="n_independent_loci >= 3 and partner_conservation >= 0.5",
                kill_reason=(
                    "the gene association is not reproduced across independent "
                    "loci, so adjacency is as likely to be coincidence"
                ),
            ),
            Gate(
                gate_id="has_testable_component",
                question=(
                    "Is there a conserved partner or a well-formed, co-located "
                    "non-coding component to propose a function around?"
                ),
                expr=(
                    "partner_conservation >= 0.5 or "
                    "(array_quality >= 0.4 and array_operon_gap_bp <= 3000)"
                ),
                kill_reason=(
                    "the locus has neither a conserved partner nor a "
                    "well-formed array within the transcriptional unit: an "
                    "unusual context on its own supports no proposal"
                ),
                rationale=(
                    "Novelty is a statement about what the locus is not. "
                    "Something has to be present for a hypothesis to be about, "
                    "and it has to be close enough to be part of the same "
                    "system -- an array ten kilobases away is a neighbour, not "
                    "a component."
                ),
            ),
            Gate(
                gate_id="not_redundant",
                question="Are the supporting sequences independent observations?",
                expr="n_clusters >= 3",
                kill_reason=(
                    "the supporting sequences collapse into fewer than three "
                    "clusters: this is one observation reported many times"
                ),
            ),
        ],
        scores=[
            Score(
                score_id="architectural_novelty",
                question="How far is this from the nearest described system?",
                expr="novelty",
                weight=2.0,
                rationale="The distance from the catalogue is the reason to look at all.",
            ),
            Score(
                score_id="partner_conservation",
                question="How consistently does the partner travel with the enzyme?",
                expr="partner_conservation",
                weight=2.0,
                rationale="Conserved adjacency across distant genomes implies selection.",
            ),
            Score(
                score_id="association_significance",
                question="How strong is the enrichment of the partner beside the enzyme?",
                expr="clamp(association_neglog_p / 6.0)",
                weight=1.5,
            ),
            Score(
                score_id="taxonomic_spread",
                question="Is the system distributed across distant lineages?",
                expr="clamp(n_phyla / 3.0) * (0.4 if clade_restricted else 1.0)",
                weight=1.0,
                rationale=(
                    "A system confined to one clade may be a recent local "
                    "accident; one spread across phyla is under selection."
                ),
            ),
            Score(
                score_id="noncoding_component",
                question="Is there a well-formed non-coding component?",
                expr="array_quality",
                weight=1.5,
                rationale=(
                    "A regularly spaced array with distinct spacers is a strong, "
                    "specific architectural signal and is directly testable."
                ),
            ),
            Score(
                score_id="testability",
                question="Can the proposal be falsified at the bench in one experiment?",
                expr="clamp(n_predictions / 3.0)",
                weight=1.0,
                rationale=(
                    "A hypothesis that cannot be killed by one experiment is not "
                    "worth the bench time, however interesting it reads."
                ),
            ),
        ],
        guidance=[
            "State the alternative explanation for every claim, and say what "
            "would distinguish it from the proposal.",
            "Prefer one candidate with a decisive experiment over five with "
            "suggestive statistics.",
        ],
    )
