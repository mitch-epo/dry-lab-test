"""The harness: many sessions in parallel, every one of them recorded.

What the harness is responsible for, beyond a thread pool:

* **isolation** -- one session's failure does not take down the batch; a failed
  session is returned as a failure and the campaign decides what that means;
* **provenance** -- every session's instructions, inputs and output land on
  disk before the result is used;
* **accounting** -- wall-clock, attempts and tokens per role, so the cost of a
  campaign is a number rather than a feeling; and
* **determinism where it is available** -- identical sessions are answered
  once and reused, which matters when a campaign re-asks the same question of
  hundreds of loci.
"""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from .backends import Backend
from .session import Session, SessionResult, TranscriptStore


@dataclass
class HarnessStats:
    """Accounting for one batch."""

    sessions: int = 0
    failures: int = 0
    cached: int = 0
    wall_s: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    by_role: dict[str, int] = field(default_factory=lambda: defaultdict(int))

    def summary(self) -> str:
        parts = [f"{self.sessions} sessions"]
        if self.cached:
            parts.append(f"{self.cached} reused")
        if self.failures:
            parts.append(f"{self.failures} failed")
        parts.append(f"{self.wall_s:.1f}s wall")
        if self.output_tokens:
            parts.append(f"{self.input_tokens}+{self.output_tokens} tokens")
        return ", ".join(parts)


class Harness:
    """Runs sessions across a pool of workers."""

    def __init__(
        self,
        backend: Backend,
        *,
        transcripts: Path | TranscriptStore | None = None,
        max_workers: int = 8,
        reuse_identical: bool = True,
    ):
        self.backend = backend
        self.max_workers = max_workers
        self.reuse_identical = reuse_identical
        if isinstance(transcripts, TranscriptStore):
            self.transcripts: TranscriptStore | None = transcripts
        elif transcripts is not None:
            self.transcripts = TranscriptStore(transcripts)
        else:
            self.transcripts = None
        self._memo: dict[str, SessionResult] = {}
        self.stats = HarnessStats()

    def run_one(self, session: Session) -> SessionResult:
        return self.run([session])[0]

    def run(
        self,
        sessions: Iterable[Session],
        *,
        on_result: Callable[[SessionResult], None] | None = None,
    ) -> list[SessionResult]:
        """Run a batch, preserving input order in the returned list."""
        import time

        sessions = list(sessions)
        results: list[SessionResult | None] = [None] * len(sessions)
        t0 = time.monotonic()

        def reuse(target: Session, cached: SessionResult) -> SessionResult:
            # Re-key the reused result to this session so the transcript and
            # any downstream lookup still resolve by session id.
            self.stats.cached += 1
            return SessionResult(
                session_id=target.session_id, role=cached.role, output=cached.output,
                backend=cached.backend + "+reused", ok=cached.ok, error=cached.error,
                attempts=cached.attempts, duration_s=0.0,
                fingerprint=cached.fingerprint,
            )

        # Deduplicate within the batch as well as against earlier batches. A
        # campaign asks the same question of many loci, and duplicates inside
        # one batch are the common case -- checking only the cross-batch memo
        # would let every one of them through, because nothing is in the memo
        # until the batch finishes.
        pending: list[tuple[int, Session]] = []
        first_of: dict[str, int] = {}
        followers: dict[int, list[int]] = {}
        for i, s in enumerate(sessions):
            if not self.reuse_identical:
                pending.append((i, s))
                continue
            if s.fingerprint in self._memo:
                results[i] = reuse(s, self._memo[s.fingerprint])
            elif s.fingerprint in first_of:
                followers.setdefault(first_of[s.fingerprint], []).append(i)
            else:
                first_of[s.fingerprint] = i
                pending.append((i, s))

        if pending:
            workers = max(1, min(self.max_workers, len(pending)))
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(self._execute, s): i for i, s in pending}
                for fut in as_completed(futures):
                    i = futures[fut]
                    try:
                        res = fut.result()
                    except Exception as exc:
                        s = sessions[i]
                        res = SessionResult(
                            session_id=s.session_id, role=s.role, output={},
                            backend=getattr(self.backend, "name", "?"), ok=False,
                            error=f"{type(exc).__name__}: {exc}",
                            fingerprint=s.fingerprint,
                        )
                    results[i] = res

        for leader, group in followers.items():
            source = results[leader]
            assert source is not None
            for i in group:
                results[i] = reuse(sessions[i], source)

        for i, s in enumerate(sessions):
            res = results[i]
            assert res is not None
            if self.transcripts is not None:
                self.transcripts.write(s, res)
            if self.reuse_identical and res.ok:
                self._memo.setdefault(s.fingerprint, res)
            self.stats.sessions += 1
            self.stats.by_role[res.role] += 1
            self.stats.input_tokens += res.input_tokens
            self.stats.output_tokens += res.output_tokens
            if not res.ok:
                self.stats.failures += 1
            if on_result is not None:
                on_result(res)

        self.stats.wall_s += time.monotonic() - t0
        return [r for r in results if r is not None]

    def _execute(self, session: Session) -> SessionResult:
        return self.backend.run(session)
