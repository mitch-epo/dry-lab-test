# Molecular design studio

A working replica of the dry-lab workflow described in Anthropic's
[*Claude discovers a novel enzyme system*](https://www.anthropic.com/news/claude-discovers-novel-enzyme-system),
built from the "How we work" section:

> Claude reads the relevant literature and reproduces the established results
> from public data to check its methods. It then searches for family members or
> genomic neighbors that fit no described system, and writes a short,
> human-readable report for each candidate that proposes a function and
> describes the evidence supporting its claims. In follow-up analyses, Claude
> critically evaluates the evidence — typically most candidates are eliminated
> at this stage. A survey may end with a single candidate worth testing, or
> with none.

That paragraph is a pipeline. This repository implements it end to end: seven
stages, a parallel session harness, and a versioned rubric that the campaign
rewrites from its own output.

```
survey ──▶ reproduce ──▶ search ──▶ report ──▶ triage ──▶ bench
           (GATE)                                  │
             ▲                                     ▼
             └───────────────── taste ◀────────────┘
                        (rewrites the rubric)
```

---

## Quick start

Python 3.10+. No API key and no network access are needed for the default
configuration.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .

studio build-corpus campaigns/rt-phage/corpus   # offline corpus (simulated)
studio run campaigns/rt-phage                   # the full campaign
studio score campaigns/rt-phage                 # grade it against withheld truth
studio show campaigns/rt-phage                  # the run summary
```

Measured on a 4-core cloud container (Python 3.11), from a fresh clone:

| step | time |
|---|---|
| `build-corpus` | 5 s |
| `run` (first time) | 5 min 30 s |
| `run` (subsequent) | 2 min 40 s |
| `pytest` | 2 min 10 s |

The first run is slower because it bootstraps a profile for each of the 24
reference families and caches them under `campaigns/rt-phage/.cache`. The work
is CPU-bound and mostly single-threaded, so clock speed matters more than core
count; `max_workers` in `campaign.yaml` parallelises agent sessions, which are
the cheap part under the default `rubric` backend.

For a faster loop while changing things, build a smaller corpus — the planted
systems and all six decoy classes are still present:

```bash
studio build-corpus campaigns/rt-phage/corpus --scale 0.35 --background 25
```

### What a run produces

```
campaigns/rt-phage/runs/<run-id>/
├── manifest.json        provenance: source, backend, rubric version, every stage
├── rubric.yaml          the rubric this run was judged under (pinned)
├── survey.json          family profile, catalogue, discriminating features
├── reproduce.txt        the methods check, and whether the gate passed
├── search.json          the funnel, the architecture groups, the rejections
├── candidates.json      candidates with every measurement attached
├── reports/*.md         one short human-readable report per candidate
├── verdicts.json        disposition, failed gates, kill reason, overclaims
├── bench/*.md           wet-lab handoff, one per candidate *system*
├── taste.txt            what separated promoted from set-aside
├── rubric.next.yaml     the next rubric version, if anything was learned
└── transcripts/*.json   every agent session, with its inputs and instructions
```

---

## The seven stages

### 1. `survey` — read the literature, build the family model

Produces two things. A **family profile**: a PSSM with Henikoff sequence
weighting, calibrated against shuffled decoys so every later hit is quoted as a
z-score against an explicit null rather than as a bare number. And a
**catalogue of described systems** (`src/studio/assets/known_systems.yaml`),
which turns prose about retrons, group II introns, DGRs, RT-Cas1 and Abi RTs
into machine-checkable architectural predicates.

The catalogue is the part that is easy to skip and expensive to skip. Without
an explicit, falsifiable statement of what is already known, "fits no described
system" has no meaning, and a survey will confidently rediscover retrons.

### 2. `reproduce` — check the methods, or stop

A **gate**, not a step. Before the pipeline may propose anything new it must
re-find what is already known, using the same code path it will use for
discovery. Three checks:

- **Positive controls.** For each described system, an exhaustive reference
  search establishes the loci it should cover, and the production path is
  scored for recall and precision against it. The two paths differ on purpose —
  permissive single-sequence searches versus the calibrated profile at the
  campaign's threshold — so agreement is informative rather than tautological.
- **Negative control on sequence.** Composition-matched shuffles must not pass
  the profile threshold. If they do, the threshold is measuring composition.
- **Negative control on the array detector.** Shuffled contig DNA must not
  yield regularly spaced arrays. If it does, every array claim downstream is an
  artefact.

Fail any of these and the run stops without proposing a candidate.

**This gate earns its keep.** During development it caught a real bug: gene
labelling used a fixed bit-score cut-off, and bit score scales with alignment
length, so a 130-residue partner family could never reach the threshold however
well it matched. Every short partner silently became an "unknown" — which reads
downstream as *novelty*. The failure mode of an insensitive labeller is not
missing candidates, it is inventing them. DGR recall sat at 0.39 and the gate
refused to let the run continue. The fix was to label with bootstrapped,
PSI-BLAST-style family profiles scored as z-scores against each profile's own
null. Recall went to 0.91 with zero false labels.

### 3. `search` — find members that fit no described system

A three-stage funnel (k-mer prefilter → ungapped extension → gapped alignment),
then genomic context, then novelty:

- each hit's neighbourhood is pulled and its genes labelled by homology, leaving
  genuinely unmatched genes as `unk`;
- each locus is scored against every catalogue entry, so novelty is a measured
  quantity with a named runner-up;
- loci are grouped by **architecture**, and per-group statistics decide whether
  an association is real: partner conservation across independent loci, Fisher
  enrichment against background loci from the same genomes with
  Benjamini–Hochberg correction, sequence redundancy, taxonomic spread.

A single locus can never support a claim about a system. A recurrent
architecture can.

### 4. `report` — one short write-up per candidate

One session per candidate, run in parallel. Each sees only its own evidence
bundle, so a report cannot borrow confidence from a sibling it never saw. Every
claim carries the measurement behind it, the method that produced it, and the
records it came from. Contradicting evidence is rendered in the same list as
supporting evidence, not in a footnote.

The report stage decides nothing. It proposes and it cites.

### 5. `triage` — critically evaluate, and eliminate

Every candidate goes through the rubric's gates, each of which names the reason
it kills. A reviewer that wants to keep a candidate must point at a passing
gate; it cannot trade a failed gate against an appealing story.

Two properties this stage is built to have:

- **It can return nothing.** A stage that always produces a winner is a ranking
  function wearing a filter's clothes.
- **It reviews the report as well as the locus.** A report can be right about
  the locus and still assert more than its evidence carries. Overclaims are
  recorded separately from gate failures, because they say something about the
  *instructions* rather than about the biology — which is what stage 7 reads.

### 6. `bench` — the wet-lab handoff

One package per candidate **system**, not per locus: fifteen loci sharing an
architecture are fifteen instances of one hypothesis, and sending the same
experiment to the bench fifteen times would inflate the apparent yield of the
survey. Each handoff leads with *what would kill the hypothesis*, then
constructs (including the active-site mutant), assays, controls and scope.

### 7. `taste` — study the hypotheses, rewrite the instructions

> Because Claude produces hypotheses so prolifically, the hypotheses themselves
> have become an object of study for us. […] What we learn goes back into the
> instructions we give Claude.

This stage mines the campaign's report corpus for features that separate
promoted from set-aside candidates (AUROC per feature), finds gates that never
fire or fire on everything, counts recurring overclaims, and emits the next
rubric version.

**The circularity problem, and what is done about it.** The verdicts were
produced *by* the rubric. Learning from them and feeding the result back into
the rubric is a loop that re-derives its own premises: features the rubric
already weights will always look discriminative, and tightening them will
always look justified. Run unexamined, this manufactures confidence.

So the stage separates findings by whether they survive that critique:

| basis | meaning | applied? |
|---|---|---|
| `external` | judged against human labels in `labels.yaml` — bench results or reviewer overrides | yes |
| `structural` | a feature the rubric does not name, an inert or saturated gate, a systematic overclaim — all measurable independently of the verdicts | yes |
| `circular` | a weight change to a criterion the rubric already uses, learned from that rubric's own verdicts | **no**, reported only |

Every proposal records its basis, and the emitted rubric records why each
change was made.

---

## The rubric is the interface

The rubric is data, not prose buried in a prompt (`src/studio/rubric.py`,
serialised to `rubric.yaml`). It has **gates** — boolean conditions whose
failure eliminates a candidate and names the reason — and **scores** —
weighted expressions in 0..1 that rank survivors.

```bash
studio rubric campaigns/rt-phage            # list gates and scores
studio rubric campaigns/rt-phage --prompt   # render as reviewer instructions
```

The same object renders to the prompt an LLM reviewer receives *and* executes
directly in the deterministic backend, so the two cannot drift apart and a
taste-stage revision reaches both at once.

Expressions are evaluated through a restricted AST walker — no `eval`, no
attribute access, no imports, a whitelist of functions. The taste stage
rewrites this file automatically, so it must never be a path to code execution.
There are tests for that.

The nine gates, each corresponding to a way a comparative-genomics hypothesis
routinely fails:

| gate | eliminates when |
|---|---|
| `catalytic_core` | the catalytic triad is disrupted — a profile hit says a protein is in a family, not that it still works |
| `domain_coverage` | the hit misses the profile's high-information columns: a fragment, not a member |
| `not_low_complexity` | the score reflects amino-acid composition rather than homology |
| `significance` | the family assignment is not significant at the campaign threshold |
| `novelty` | the locus is accounted for by a described system |
| `partner_is_not_known` | the "unknown" partner is a diverged member of a described family |
| `has_testable_component` | neither a conserved partner nor a well-formed, co-located array — an unusual context alone supports no proposal |
| `reproducible_association` | the association is not reproduced across independent loci |
| `not_redundant` | the supporting sequences collapse into one observation reported many times |

---

## The agent harness

> …sometimes with a harness of our own that coordinates many Claude sessions
> running in parallel.

`studio.agents` runs sessions across a thread pool with isolation (one failure
does not take down the batch), provenance (instructions, inputs and outputs on
disk before the result is used), accounting, and deduplication of identical
work both within and across batches.

Three backends behind one interface:

- **`rubric`** (default) — executes the rubric in code. Deterministic, free,
  always available. Not a stand-in for judgement; it *is* the rubric, run. Since
  the rubric is also what an LLM reviewer is handed, this backend is the
  reference implementation of the standard the model is asked to meet, and
  running a campaign both ways measures the difference.
- **`anthropic[:model]`** — real Claude sessions through the Messages API, with
  the rubric as the system prompt and output pinned by a tool schema. Needs
  `pip install -e '.[live]'` and `ANTHROPIC_API_KEY`.
- **`replay`** — replays recorded transcripts by fingerprint, refusing any whose
  inputs have changed, so a published campaign can be re-derived exactly.

```bash
studio run campaigns/rt-phage --backend anthropic:claude-opus-5
studio run campaigns/rt-phage --backend replay
```

---

## Data sources

Everything reads through one interface (`studio.sources.base.DataSource`), and
the run manifest records which was used.

**`live`** — UniProtKB, NCBI E-utilities, InterPro/Pfam and Europe PMC, all
public and unauthenticated, with rate limiting and pagination.

**`local:<path>`** — a generated corpus, for running where those hosts are not
reachable.

> **The environment this was developed in blocks all four hosts at the egress
> proxy**, so the live backend is written against the documented REST contracts
> but has **not been executed end to end**. `studio doctor` probes reachability
> and reports which backends are usable, so a live run fails at startup rather
> than silently halfway through.

### The offline corpus is simulated, and says so everywhere

Every sequence in `campaigns/rt-phage/corpus` is generated by the evolution
simulator in `studio/corpus/evolve.py`. **None is a natural sequence and no
biological conclusion may be drawn from it.** The provenance block is marked
synthetic, and the banner is propagated into the manifest and onto the first
line of every rendered report.

For the corpus to be worth anything as a test, the homology in it has to be
*earned*. Families are evolved along random trees under a substitution model
whose exchangeabilities come from BLOSUM62, with site-specific rates: catalytic
motifs invariant, structural cores slow, loops fast with indels. The result has
the properties real families have — graded identity (RT members sit at 28–42%
pairwise), conserved motifs inside variable backgrounds — and is therefore a
fair test of a homology search rather than a rehearsal of the generator.

---

## Scoring: does the pipeline actually work?

`studio score` grades a run against ground truth the corpus withholds from
every stage (`TRUTH.withheld.json`, loaded only by the scoring command). Three
things are planted, and the pipeline is told about only the first:

1. instances of the **described** systems — the reproduce stage must recover
   these;
2. two genuinely **novel** architectures, one of them mirroring the article's
   system (enzyme + partner gene of unknown function + a long array of evenly
   spaced repeats, in phage genomes);
3. six classes of **decoy**, each built to be attractive to a naive search and
   killed by one specific check.

The decoys are the point. A pipeline that finds the planted systems but also
promotes the decoys has demonstrated nothing, because this workflow is defined
by what it discards.

Current result on the reference corpus:

```
Novel systems planted   : 25
  reached the candidate list: 23
  promoted to the bench     : 23   (resolving to 2 candidate systems)

[CAUGHT] decoy:dead-rt            killed by catalytic_core
[CAUGHT] decoy:rt-fragment        never reached triage (core coverage, at search)
[CAUGHT] decoy:distant-array      killed by has_testable_component
[CAUGHT] decoy:vntr               killed by has_testable_component
[CAUGHT] decoy:diverged-cas1      killed by partner_is_not_known
[CAUGHT] decoy:singleton-partner  killed by reproducible_association

Described systems promoted: 0
VERDICT: clean
```

Every decoy is killed by the check it was designed to test, and the top-ranked
system is the ART-like architecture. Note that a survival rate this high
reflects a corpus deliberately enriched with planted systems; on real data the
article's expectation — most candidates eliminated, a survey ending with one
candidate or none — is the one to expect.

Two decoys are worth describing because an earlier version of them was wrong.
`decoy:vntr` and `decoy:distant-array` originally carried the *same* conserved
partner family as the planted novel system, which made them genuine
RT-plus-conserved-unknown loci that merely happened to have a bad array —
promoting those was arguably correct, so they tested nothing. A decoy is only a
decoy when the feature under test is its sole apparent evidence.

---

## Tests

```bash
pip install -e '.[dev]'
pytest          # 68 tests, about two minutes
```

The suite pins the properties conclusions rest on: BLOSUM62 values and
alignment correctness; that batched scanning agrees with the pairwise aligner
*exactly* (the optimisation must change speed and nothing else); that motif and
array detection distinguish the cases triage is asked to distinguish; that
Fisher and Benjamini–Hochberg produce known values; that the rubric evaluator
refuses arbitrary code; that the harness isolates failures; and an end-to-end
campaign asserting that the planted system is recovered and every decoy class
is eliminated.

---

## Notable implementation details

**Anti-diagonal Gotoh.** On an anti-diagonal every cell depends only on the two
preceding ones, so the whole affine-gap recurrence — including the horizontal
gap state, the part that blocks naive row-wise vectorisation — becomes three
numpy operations per diagonal.

**Batched scanning (`profile_scores_batch`).** The single-pair DP spends almost
all its time in interpreter overhead, not arithmetic. Batching runs every
target through the same anti-diagonal together, so the loop runs once for the
whole batch: 8× faster with bit-identical scores, and full alignments on the
reference corpus dropped from 212 to 15 because scoring now precedes traceback.

**Array-class grouping.** An early version keyed architecture groups on exact
array geometry (`array[21x23]`). Copy number and repeat length vary between
instances of one real system, so every ART-like locus formed a group of one and
was then eliminated for non-reproducibility — the component that made the
system most recognisable was what disqualified it. Grouping is now on
architectural *class* (`array[spacer-array]`), with geometry kept for the
report.

**Motif placement as constraint satisfaction.** Choosing each motif's first
match independently lets a permissive upstream pattern matching by chance block
the genuine downstream match of an essential one, and the set then reads as
disrupted when it is intact. Placement searches for the best mutually
consistent assignment, ranked by essential motifs placed, then total placed,
then pattern information content.

---

## Layout

```
src/studio/
├── bio/           analysis kernel — alignment, profiles, motifs, arrays,
│                  genomic context, clustering, statistics. Network-free
│                  and deterministic.
├── sources/       DataSource interface; live (UniProt/NCBI/InterPro/EuropePMC)
│                  and local backends
├── corpus/        the evolution simulator and the corpus builder
├── agents/        session, harness, backends (rubric / anthropic / replay)
├── stages/        the seven stages
├── assets/        described-systems catalogue, literature digest
├── rubric.py      gates, scores, the restricted expression evaluator
├── schemas.py     Evidence, Candidate, CandidateReport, Verdict, Manifest
├── pipeline.py    stage sequencing and the reproduce gate
├── scoring.py     grading against withheld ground truth
└── cli.py         the studio command
```

## Commands

| command | does |
|---|---|
| `studio build-corpus <dir>` | generate the offline corpus |
| `studio run <campaign>` | run a campaign (`--backend`, `--stages`, `--force`) |
| `studio score <campaign>` | grade against withheld ground truth |
| `studio show <campaign>` | print the summary, or `--candidate` for one report |
| `studio rubric <campaign>` | show the rubric, or `--prompt` to render it |
| `studio list-runs <campaign>` | all runs, with gate status and counts |
| `studio doctor` | what this environment can reach, and which backends work |
