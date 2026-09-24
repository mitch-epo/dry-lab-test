> SIMULATED DATA -- every sequence in this run is generated. No biological claim in these outputs refers to a real organism.

# RT with a conserved partner and a regularly spaced non-coding array

**cand-ctg00191_g001**

**Triage: PROMOTE** (priority 0.88, rubric v1)

## Proposed function

A three-component retroelement: the reverse transcriptase, a conserved partner protein of unknown function encoded in the same operon, and an adjacent array of evenly spaced repeats that is part of the same locus. The architecture matches no described RT-associated system. The most economical reading is that the array is the RT's template or substrate -- as msr/msd is for a retron and TR is for a DGR -- and that the partner is required for the reaction or is its effector, in a phage-encoded defence or counter-defence context. Observed architecture: RT(+)-unk(+)-array[21x23].

## Summary

cand-ctg00191_g001: reverse transcriptase ctg00191_g001 on ctg00191 (unclassified-siphovirus). Novelty 1.00 against the catalogue. The nearest described system on gene content alone is Group II intron RT (maturase) (1.00), but that match is ruled out on architecture -- which is why the novelty score is high despite the gene-content similarity. 5 supporting and 1 contradicting observations. The architecture recurs at 20 loci (20 non-redundant, 1 phyla).

## Evidence

- [direct, supports] catalytic triad intact: D117 (motif A LVDI@115), D195/D196 (motif C YADD@193)
  - method: PROSITE-style motif model of the RT palm domain (motifs A and C)
  - records: `ctg00191_g001`
- [direct, supports] profile hit at 92 bits (E=4.3e-20, z=113 above the shuffled-decoy null), covering 88% of the profile's core columns
  - method: calibrated PSSM of 334 columns vs. the protein database
  - records: `ctg00191_g001`
- [strong, supports] locus architecture is RT(+)-unk(+)-array[21x23]; the closest described system is Group II intron RT (maturase), which the gene content matches at 1.00, but it is ruled out by unexpected repeat array
  - method: operon prediction plus architectural comparison to the described-system catalogue
  - records: `ctg00191`
- [strong, supports] the partner travels with the enzyme at 100% of the 20 loci sharing this architecture; OR=266.2, p=2.39e-16 (with/without: 20/0 in set, 33/217 in background)
  - method: Fisher exact test against background loci from the same genomes, BH-corrected
  - records: `ctg00191_g002`
- [strong, supports] 21x 23 bp repeat, spacers 44 bp (regularity 0.97), repeat identity 0.98, spacer identity 0.25, span 4396-5750 -- CRISPR-like (repeat-spacer array with distinct spacers)
  - method: repeat-spacer array detection on the contig, with spacer-distinctness scoring
  - records: `ctg00191`
- [weak, **against**] 20 sequences -> 20 clusters at 90% id across 3 taxa / 1 phyla; largest clade holds 100%
  - method: greedy clustering at 90% identity plus lineage tabulation
  - records: `ctg00206_g001`, `ctg00183_g002`, `ctg00192_g001`, `ctg00191_g001`, `ctg00181_g001`

## Alternative explanations

- The locus is a diverged Group II intron RT (maturase) (gene content matches at 1.00) whose partner has fallen below the labelling threshold.
- The array is an unrelated element that happens to lie nearby; the nearest one is 168 bp from the predicted operon.
- The architecture is a recent, clade-local arrangement rather than a conserved system.
- The RT is a mobile element that has inserted next to an unrelated gene and array, with no functional relationship between them.

## Predictions and how to falsify them

1. The purified protein has RNA-dependent DNA polymerase activity requiring the motif A and motif C aspartates.
   - assay: Express and purify the protein; assay primer extension on an RNA template with Mg2+, alongside a motif C D->N active-site mutant.
   - would falsify: No extension by the wild type, or equal extension by the active-site mutant.
2. The array is transcribed and processed into discrete short RNAs with boundaries at the repeats.
   - assay: Small-RNA sequencing and Northern blot on the native host, or on a heterologous host carrying the cloned locus.
   - would falsify: No transcript over the array, or a single unprocessed transcript with no repeat-aligned ends.
3. The array is the RT's template: cDNA products map to array-derived sequence.
   - assay: Sequence nucleic acid co-purifying with the RT; map reads to the locus and check strand and endpoints against the repeat units.
   - would falsify: Co-purifying nucleic acid maps elsewhere, or no DNA product is recoverable.
4. The partner protein is required for the reaction, or forms a complex with the RT.
   - assay: Co-express and pull down; repeat the polymerase assay with and without the partner.
   - would falsify: No co-purification and no change in activity when the partner is omitted.
5. The locus affects phage-host outcome: expressing it changes plating efficiency on challenge.
   - assay: Clone the locus into a permissive host; measure efficiency of plaquing against a panel of phages, with an empty-vector control and an active-site mutant.
   - would falsify: No change in efficiency of plaquing for any phage tested.

## Open questions

- Is the array transcribed from one promoter or many?
- Do the spacers share a source, or are they locus-specific?
- Does the partner have a fold that suggests a biochemical role?
- Is the narrow taxonomic distribution real, or an artefact of which genomes have been sequenced?

---
_Confidence 0.77; written by rubric-reporter; rubric v1; 2026-09-24T09:38:33+00:00._