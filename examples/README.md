# Example output

Rendered artifacts from one reference run of the `rt-phage` campaign, committed
so the repository shows what the pipeline produces without requiring a run.

Reproduce them with:

```bash
studio build-corpus campaigns/rt-phage/corpus
studio run campaigns/rt-phage --run-id demo
studio score campaigns/rt-phage --run-id demo
```

The corpus is generated from a fixed seed, so a rebuild gives the same data and
the same result.

| file | what it is |
|---|---|
| `summary.md` | the run summary: every stage, what it did, and the outcome |
| `reproduce.txt` | the methods check — positive and negative controls, and the gate |
| `report-promoted.md` | a candidate report that survived triage |
| `report-eliminated.md` | one that did not, with the check that killed it |
| `bench-handoff.md` | the wet-lab package for the top candidate system |
| `taste.txt` | what separated promoted from set-aside, and the rubric changes |
| `score.txt` | grading against the corpus's withheld ground truth |

Everything here is **simulated data** — see the banner on each report. The
corpus exists so the pipeline can be run and scored offline; no biological
claim in these files refers to a real organism.
