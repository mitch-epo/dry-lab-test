# Working in this repository

A dry-lab studio that turns an uncharacterised protein family into a short list
of lab-testable hypotheses. See README.md for what it does and why.

## Ground rules

**The corpus is simulated.** Every sequence under `campaigns/*/corpus` is
generated. Never present a result from it as a biological finding, and never
remove the `SIMULATED DATA` banner from a rendered report or manifest.

**Never read `TRUTH.withheld.json` from a pipeline stage.** It is the answer
key. Only `studio score` may open it. A stage that consults it silently
invalidates every number the project produces.

**Do not tune thresholds to make the reproduce gate pass.** The gate exists to
catch methodological errors, and it has caught real ones. When it fails, find
the cause. Two fixes so far came from doing that: labelling by length-aware
z-scores instead of raw bit score, and grouping loci by architecture class
instead of array geometry. Both were bugs the gate surfaced; neither would have
been found by loosening a number.

**A change that makes the pipeline promote more is not automatically an
improvement.** Run `studio score` after any change to the search, rubric or
triage logic. Promoting a decoy is a regression even if the planted systems are
also found.

## Verifying a change

```bash
pytest                                    # 68 tests, ~3 min
studio run campaigns/rt-phage --run-id check
studio score campaigns/rt-phage --run-id check    # must print "VERDICT: clean"
```

`studio score` exits non-zero when a decoy leaks or a described system is
promoted, so it works as a check in CI.

The labeller is cached under `campaigns/*/.cache`; delete it after changing
anything in `bio/profile.py`, `bio/iterate.py` or the reference set.

## Conventions

**Do not use `X | Y` unions inside a pydantic model.** Pydantic resolves a
model's annotations at runtime, so PEP 604 syntax is evaluated even under
`from __future__ import annotations` and raises on Python 3.9 — which is the
floor, because it is what stock macOS ships. Use `Optional[X]` / `Union[X, Y]`
in `schemas.py`, and anywhere else a library reads annotations at runtime
(pydantic models, typer command signatures). Ordinary function signatures use
`|` freely; those stay strings and are never evaluated. The failure is silent
on a newer interpreter, so CI runs the suite on 3.9 as well as on the newest
release.

- `bio/` is network-free and deterministic. Keep it that way: the reproduce
  stage is meaningless if methods drift between runs.
- Every claim a report makes carries an `Evidence` record with a `method` and
  the records it came from. Do not add prose claims without one.
- Rubric expressions are evaluated by a restricted AST walker, never `eval`.
  Do not widen it.
- New gates state what they eliminate and why, in the gate's own
  `kill_reason` and `rationale`. A gate whose failures cannot be explained to
  a reader is not usable.
