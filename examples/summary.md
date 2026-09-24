# Campaign run demo

> SIMULATED DATA -- every sequence in this run is generated. No biological claim in these outputs refers to a real organism.

- campaign: **rt-phage**
- source: local (generated corpus at /home/user/dry-lab-test/campaigns/rt-phage/corpus)
- backend: rubric
- rubric: v1
- reproduce gate: PASSED

## Stages

### survey (1.7s)
- catalogue v1 with 6 described systems
- profile built from 14 seed sequences, 102 core columns above 1.5 bits
- decoy calibration: mu=15.9 sigma=1.9 from shuffled seeds
- loaded 24 family profiles for gene labelling

### reproduce (144.9s)
- production path recovered 219 anchor loci at z>=8.0, E<=1e-05
- reference path recovered 235 anchor loci at E<=0.001
- labels: strict z>=6.0, sensitive z>=4.0
- 4 positive controls, 2 negative controls

### search (67.8s)
- 1890 in database -> 1869 past k-mer prefilter -> 280 past ungapped extension -> 230 aligned -> 219 reported
- 219 loci assembled; 72 scored novelty >= 0.35
- 8 distinct architectures among the novel loci
- 250 background loci sampled for enrichment tests

### report (0.1s)
- 60 sessions, 0.1s wall

### triage (0.1s)
- 120 sessions, 0.2s wall
- promote 23, hold 0, eliminate 37 of 60
- most frequent gate failures: partner_is_not_known x20, reproducible_association x13, has_testable_component x9, catalytic_core x6, not_redundant x4

### bench (0.0s)
- 2 candidate system(s) from 23 promoted loci (2 distinct architectures)

### taste (0.0s)
- 60 reports analysed; labels from verdicts
- 2 of 7 proposals applied
- labels came from the same rubric being revised, so weight changes to existing criteria are reported but NOT applied; supply labels.yaml with external outcomes to lift that restriction

## Outcome

60 candidates written up, 23 promoted, resolving to 2 candidate system(s) for the bench.

### cand-ctg00191_g001 -- RT with a conserved partner and a regularly spaced non-coding array

- architecture: `RT(+)-unk(+)-array[21x23]`
- independent loci: 14
- report: `reports/cand-ctg00191_g001.md`
- bench handoff: `bench/cand-ctg00191_g001.md`

### cand-ctg00196_g002 -- RT with a conserved partner of unknown function

- architecture: `RT(+)-unk(+)`
- independent loci: 9
- report: `reports/cand-ctg00196_g002.md`
- bench handoff: `bench/cand-ctg00196_g002.md`

## Taste analysis

```
Taste analysis
==============

corpus: 60 reports, 23 promoted, 37 set aside
labels: verdicts

Feature discrimination (ranked by distance from chance):
  partner_known_neglog_e: AUROC 0.24 (in rubric); promoted mean 0.0291 vs set-aside 14.9
  partner_known_homology_bits: AUROC 0.24 (NOT in rubric); promoted mean 1.83 vs set-aside 61
  array_spacer_identity: AUROC 0.25 (NOT in rubric); promoted mean 0.145 vs set-aside 0.325
  novelty: AUROC 0.32 (in rubric); promoted mean 0.535 vs set-aside 0.732
  n_phyla: AUROC 0.34 (in rubric); promoted mean 1.78 vs set-aside 3
  n_predictions: AUROC 0.65 (in rubric); promoted mean 3.83 vs set-aside 3.32
  partner_conservation: AUROC 0.62 (in rubric); promoted mean 1 vs set-aside 0.757
  association_significant: AUROC 0.62 (NOT in rubric); promoted mean 1 vs set-aside 0.757
  n_independent_loci: AUROC 0.60 (in rubric); promoted mean 15.7 vs set-aside 12
  n_clusters: AUROC 0.60 (in rubric); promoted mean 15.7 vs set-aside 12
  triad_intact: AUROC 0.58 (in rubric); promoted mean 1 vs set-aside 0.838
  clade_restricted: AUROC 0.57 (in rubric); promoted mean 0.609 vs set-aside 0.46

Gate activity:
  partner_is_not_known: eliminated 20
  reproducible_association: eliminated 13
  has_testable_component: eliminated 9
  catalytic_core: eliminated 6
  not_redundant: eliminated 4
  inert (never fired): domain_coverage, not_low_complexity, significance, novelty

Proposals:
  [held/circular] reweight partner_known_neglog_e: separation AUROC 0.24 suggests raising this criterion's weight
  [APPLIED/structural] add_score partner_known_homology_bits: add as a ranking criterion with weight 0.5 (inverted: lower is better)
  [APPLIED/structural] add_score array_spacer_identity: add as a ranking criterion with weight 0.5 (inverted: lower is better)
  [held/structural] retire_gate domain_coverage: eliminated nothing in this corpus; either it is redundant with another gate or its threshold is unreachable
  [held/structural] retire_gate not_low_complexity: eliminated nothing in this corpus; either it is redundant with another gate or its threshold is unreachable
  [held/structural] retire_gate significance: eliminated nothing in this corpus; either it is redundant with another gate or its threshold is unreachable
  [held/structural] retire_gate novelty: eliminated nothing in this corpus; either it is redundant with another gate or its threshold is unreachable
```
