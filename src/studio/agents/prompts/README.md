# Prompts

There are no prompt template files here, and that is deliberate.

The instructions an agent session receives are assembled in code from two
pieces:

* the role's standing instructions, held next to the stage that issues them
  (`stages/report.py` and `stages/triage.py` each define one); and
* `Rubric.as_prompt()`, which renders the campaign's current gates, scores,
  thresholds and learned guidance.

Keeping the rubric as the single source means the deterministic backend and a
model backend are held to the same standard: one executes the expressions, the
other reads them, and both come from the same object. A file of frozen prompt
text beside a live rubric would drift from it within one taste cycle, and the
drift would be invisible -- the model would be judged against yesterday's
standard while the scores were computed against today's.

To see what a reviewer session actually receives:

    studio rubric campaigns/rt-phage --prompt
