"""A session: one unit of agent work, with its transcript kept.

The article describes "a harness of our own that coordinates many Claude
sessions running in parallel". A session here is the coordinated unit: a role,
an evidence bundle, an expected output shape, and a record of exactly what was
asked and what came back.

Transcripts are not a debugging convenience. A campaign that produces hundreds
of reports is only auditable if each one can be traced to the inputs and the
instruction version that produced it, which is also the raw material the taste
stage mines.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _stable_hash(obj: Any) -> str:
    blob = json.dumps(obj, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


@dataclass
class Session:
    """A request to an agent backend."""

    session_id: str
    role: str                                   # "reporter", "reviewer", "analyst"
    instructions: str                           # the system/rubric text
    task: str                                   # what to do with this bundle
    bundle: dict[str, Any] = field(default_factory=dict)   # the evidence
    schema_hint: str = ""                       # shape of the expected JSON
    max_tokens: int = 4096
    temperature: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def fingerprint(self) -> str:
        """Stable identity of this request, for caching and replay."""
        return _stable_hash(
            {"role": self.role, "instructions": self.instructions,
             "task": self.task, "bundle": self.bundle}
        )

    def render_prompt(self) -> str:
        """The user-turn text a model backend receives."""
        parts = [self.task, ""]
        if self.schema_hint:
            parts += ["Return a single JSON object with this shape:",
                      self.schema_hint, ""]
        parts += ["Evidence bundle:", json.dumps(self.bundle, indent=2, default=str)]
        return "\n".join(parts)


@dataclass
class SessionResult:
    """What a backend returned, plus the accounting."""

    session_id: str
    role: str
    output: dict[str, Any]
    backend: str
    ok: bool = True
    error: str = ""
    attempts: int = 1
    duration_s: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    raw: str = ""
    fingerprint: str = ""
    finished_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )

    def to_transcript(self, session: Session) -> dict[str, Any]:
        return {
            "session_id": self.session_id,
            "role": self.role,
            "backend": self.backend,
            "fingerprint": self.fingerprint or session.fingerprint,
            "ok": self.ok,
            "error": self.error,
            "attempts": self.attempts,
            "duration_s": round(self.duration_s, 3),
            "tokens": {"input": self.input_tokens, "output": self.output_tokens},
            "finished_at": self.finished_at,
            "instructions": session.instructions,
            "task": session.task,
            "bundle": session.bundle,
            "output": self.output,
            "raw": self.raw,
        }


class TranscriptStore:
    """Append-only transcript storage, one JSON file per session."""

    def __init__(self, directory: Path):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)

    def write(self, session: Session, result: SessionResult) -> Path:
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in session.session_id)
        path = self.directory / f"{safe}.json"
        path.write_text(json.dumps(result.to_transcript(session), indent=1, default=str))
        return path

    def read(self, session_id: str) -> dict[str, Any] | None:
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in session_id)
        path = self.directory / f"{safe}.json"
        if not path.exists():
            return None
        return json.loads(path.read_text())

    def all(self) -> list[dict[str, Any]]:
        return [json.loads(p.read_text()) for p in sorted(self.directory.glob("*.json"))]
