"""Agent backends: who actually answers a session.

Three of them, all satisfying one interface:

``rubric``
    Applies the campaign rubric to the evidence bundle in code. Deterministic,
    free, and always available. It is not a stand-in for judgement -- it *is*
    the rubric, executed. Because the rubric is also what an LLM reviewer is
    handed, this backend is the reference implementation of the standard the
    model is being asked to meet, and running a campaign both ways measures the
    difference.

``anthropic``
    Real Claude sessions through the Messages API, with the rubric as the
    system prompt and a JSON schema enforced through tool use.

``replay``
    Replays recorded transcripts by fingerprint, so a published campaign can be
    re-derived exactly without spending tokens or depending on model version.
"""

from __future__ import annotations

import json
import os
import re
import time
from typing import Any, Protocol

from ..evidence import overclaims
from ..rubric import Rubric
from .session import Session, SessionResult, TranscriptStore


class Backend(Protocol):
    name: str

    def run(self, session: Session) -> SessionResult: ...


# --- deterministic rubric backend --------------------------------------------


class RubricBackend:
    """Executes the campaign rubric against the evidence bundle."""

    name = "rubric"

    def __init__(self, rubric: Rubric):
        self.rubric = rubric

    def run(self, session: Session) -> SessionResult:
        t0 = time.monotonic()
        try:
            handler = {
                "reviewer": self._review,
                "reporter": self._report,
                "analyst": self._analyse,
            }.get(session.role)
            if handler is None:
                raise ValueError(f"no rubric handler for role {session.role!r}")
            output = handler(session)
            ok, err = True, ""
        except Exception as exc:
            output, ok, err = {}, False, f"{type(exc).__name__}: {exc}"
        return SessionResult(
            session_id=session.session_id,
            role=session.role,
            output=output,
            backend=self.name,
            ok=ok,
            error=err,
            duration_s=time.monotonic() - t0,
            fingerprint=session.fingerprint,
        )

    # -- roles --------------------------------------------------------------

    def _review(self, session: Session) -> dict[str, Any]:
        """Apply gates and scores; produce a disposition with reasons."""
        metrics = dict(session.bundle.get("metrics", {}))
        failed, reasons, scores, priority = self.rubric.apply(metrics)
        disposition = self.rubric.disposition(failed, priority)

        critique: list[str] = list(reasons)

        # Check the report's own claims against the measurements behind them.
        # This is the part of triage that is about the write-up rather than the
        # locus: a report may be right about the locus and still assert more
        # than its evidence carries. The rules come from studio.evidence, the
        # same module the reporter assigns strengths with, so a flagged
        # overclaim means the report reached past its evidence rather than that
        # the two stages disagree about the scale.
        overclaim_list = overclaims(session.bundle.get("claims", []), metrics)

        if not failed and priority < self.rubric.promote_threshold:
            critique.append(
                f"priority {priority:.2f} is below the promote threshold "
                f"{self.rubric.promote_threshold:.2f}: "
                + ", ".join(
                    f"{k} {v:.2f}" for k, v in sorted(scores.items(), key=lambda kv: kv[1])[:3]
                )
            )

        return {
            "candidate_id": session.bundle.get("candidate_id", ""),
            "disposition": disposition,
            "priority": round(priority, 4),
            "failed_gates": failed,
            "kill_reason": "; ".join(reasons),
            "scores": {k: round(v, 4) for k, v in scores.items()},
            "critique": critique,
            "overclaims": overclaim_list,
        }

    def _report(self, session: Session) -> dict[str, Any]:
        """Write a candidate report from the evidence bundle.

        Deterministic prose: every sentence is generated from a measurement, so
        the report cannot say anything the bundle does not contain. That is a
        narrower writer than a model, and a strictly honest one.
        """
        from ..render import compose_report

        return compose_report(session.bundle, self.rubric)

    def _analyse(self, session: Session) -> dict[str, Any]:
        """Meta-analysis of the report corpus -- implemented by the taste stage."""
        from ..stages.taste import analyse_corpus

        return analyse_corpus(session.bundle, self.rubric)


# --- Anthropic backend --------------------------------------------------------

_JSON_BLOCK = re.compile(r"```(?:json)?\s*(\{.*?\})\s*```", re.S)


class AnthropicBackend:
    """Real Claude sessions through the Messages API.

    The rubric goes in the system prompt; the evidence bundle goes in the user
    turn; the output shape is pinned with a tool definition so the reply is
    parsed rather than scraped. Needs ``pip install molecular-design-studio[live]``
    and ``ANTHROPIC_API_KEY``.
    """

    name = "anthropic"

    def __init__(
        self,
        model: str = "claude-opus-5",
        *,
        api_key: str | None = None,
        max_retries: int = 3,
        timeout: float = 120.0,
    ):
        try:
            import anthropic
        except ImportError as exc:  # pragma: no cover - depends on optional extra
            raise RuntimeError(
                "the anthropic backend needs the 'live' extra: "
                "pip install 'molecular-design-studio[live]'"
            ) from exc
        key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not key:
            raise RuntimeError("ANTHROPIC_API_KEY is not set")
        self.model = model
        self.max_retries = max_retries
        self._client = anthropic.Anthropic(api_key=key, timeout=timeout)

    def run(self, session: Session) -> SessionResult:
        t0 = time.monotonic()
        last_error = ""
        for attempt in range(1, self.max_retries + 1):
            try:
                msg = self._client.messages.create(
                    model=self.model,
                    max_tokens=session.max_tokens,
                    temperature=session.temperature,
                    system=session.instructions,
                    messages=[{"role": "user", "content": session.render_prompt()}],
                    tools=[{
                        "name": "submit",
                        "description": "Submit the structured result.",
                        "input_schema": {"type": "object", "additionalProperties": True},
                    }],
                    tool_choice={"type": "tool", "name": "submit"},
                )
                output: dict[str, Any] = {}
                raw_parts: list[str] = []
                for block in msg.content:
                    if getattr(block, "type", "") == "tool_use":
                        output = dict(block.input)
                    elif getattr(block, "type", "") == "text":
                        raw_parts.append(block.text)
                raw = "\n".join(raw_parts)
                if not output and raw:
                    output = self._salvage(raw)
                if not output:
                    raise ValueError("model returned no structured output")
                return SessionResult(
                    session_id=session.session_id,
                    role=session.role,
                    output=output,
                    backend=f"{self.name}:{self.model}",
                    attempts=attempt,
                    duration_s=time.monotonic() - t0,
                    input_tokens=getattr(msg.usage, "input_tokens", 0),
                    output_tokens=getattr(msg.usage, "output_tokens", 0),
                    raw=raw,
                    fingerprint=session.fingerprint,
                )
            except Exception as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                if attempt < self.max_retries:
                    time.sleep(min(2 ** attempt, 16))

        return SessionResult(
            session_id=session.session_id, role=session.role, output={},
            backend=f"{self.name}:{self.model}", ok=False, error=last_error,
            attempts=self.max_retries, duration_s=time.monotonic() - t0,
            fingerprint=session.fingerprint,
        )

    @staticmethod
    def _salvage(text: str) -> dict[str, Any]:
        """Last-resort parse of a JSON object out of prose."""
        m = _JSON_BLOCK.search(text)
        candidate = m.group(1) if m else None
        if candidate is None:
            start, depth = text.find("{"), 0
            if start >= 0:
                for i in range(start, len(text)):
                    depth += (text[i] == "{") - (text[i] == "}")
                    if depth == 0:
                        candidate = text[start:i + 1]
                        break
        if not candidate:
            return {}
        try:
            parsed = json.loads(candidate)
        except json.JSONDecodeError:
            return {}
        return parsed if isinstance(parsed, dict) else {}


# --- replay backend -----------------------------------------------------------


class ReplayBackend:
    """Replays recorded transcripts, matching on session id then fingerprint."""

    name = "replay"

    def __init__(self, store: TranscriptStore, *, strict: bool = True):
        self.store = store
        self.strict = strict
        self._by_fingerprint = {
            t.get("fingerprint", ""): t for t in store.all() if t.get("fingerprint")
        }

    def run(self, session: Session) -> SessionResult:
        record = self.store.read(session.session_id)
        if record is None:
            record = self._by_fingerprint.get(session.fingerprint)
        if record is None:
            if self.strict:
                raise KeyError(
                    f"no recorded transcript for {session.session_id} "
                    f"(fingerprint {session.fingerprint})"
                )
            return SessionResult(
                session_id=session.session_id, role=session.role, output={},
                backend=self.name, ok=False, error="no recorded transcript",
                fingerprint=session.fingerprint,
            )
        if self.strict and record.get("fingerprint") != session.fingerprint:
            raise ValueError(
                f"transcript for {session.session_id} was recorded against "
                f"different inputs (fingerprint {record.get('fingerprint')} "
                f"!= {session.fingerprint}); re-run without --backend replay"
            )
        return SessionResult(
            session_id=session.session_id,
            role=session.role,
            output=record.get("output", {}),
            backend=f"{self.name}<{record.get('backend', '?')}>",
            ok=record.get("ok", True),
            error=record.get("error", ""),
            fingerprint=record.get("fingerprint", ""),
        )


def make_backend(spec: str, rubric: Rubric, *, transcripts: TranscriptStore | None = None) -> Backend:
    """Resolve a backend specification.

    ``rubric`` | ``anthropic`` | ``anthropic:<model>`` | ``replay``
    """
    if spec == "rubric":
        return RubricBackend(rubric)
    if spec == "replay":
        if transcripts is None:
            raise ValueError("the replay backend needs a transcript store")
        return ReplayBackend(transcripts)
    if spec == "anthropic" or spec.startswith("anthropic:"):
        model = spec.split(":", 1)[1] if ":" in spec else "claude-opus-5"
        return AnthropicBackend(model)
    raise ValueError(f"unknown backend {spec!r}")
