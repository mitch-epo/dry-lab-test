> SIMULATED DATA -- every sequence in this run is generated. No biological claim in these outputs refers to a real organism.

# Bench handoff -- cand-ctg00191_g001

**RT with a conserved partner and a regularly spaced non-coding array**

Architecture `RT(+)-unk(+)-array[21x23]`, seen at 14 independent loci. This handoff covers the system; the other instances are listed at the end and are not separate experiments.

## Proposed function

A three-component retroelement: the reverse transcriptase, a conserved partner protein of unknown function encoded in the same operon, and an adjacent array of evenly spaced repeats that is part of the same locus. The architecture matches no described RT-associated system. The most economical reading is that the array is the RT's template or substrate -- as msr/msd is for a retron and TR is for a DGR -- and that the partner is required for the reaction or is its effector, in a phage-encoded defence or counter-defence context. Observed architecture: RT(+)-unk(+)-array[21x23].

## What would kill this hypothesis

- The purified protein shows no primer extension on any template under any tested condition: the family assignment is wrong or the protein is not an active enzyme.
- The active-site mutant retains the activity: the readout is a contaminant, not the enzyme.
- No transcript is detectable over the array, or the transcript is a single unprocessed species: the array is not the non-coding component the proposal requires.
- Nucleic acid co-purifying with the enzyme does not map to the array: the array is adjacent but not the template.
- The partner neither co-purifies nor changes activity: the conserved adjacency does not reflect a physical or functional partnership.
- No extension by the wild type, or equal extension by the active-site mutant.
- No transcript over the array, or a single unprocessed transcript with no repeat-aligned ends.
- Co-purifying nucleic acid maps elsewhere, or no DNA product is recoverable.
- No co-purification and no change in activity when the partner is omitted.
- No change in efficiency of plaquing for any phage tested.

## Constructs

### pEXP-anchor-WT
- insert: ctg00191_g001 coding sequence, codon-optimised
- tag: N-terminal His6-SUMO, cleavable
- host: E. coli BL21(DE3)
- purpose: purify the enzyme for in vitro polymerase assays

### pEXP-anchor-DD/NN
- insert: same, with the motif C aspartate pair mutated to asparagine
- tag: N-terminal His6-SUMO, cleavable
- host: E. coli BL21(DE3)
- purpose: active-site control: every activity attributed to the enzyme must disappear here, or it is not the enzyme's activity

### pEXP-partner
- insert: the conserved partner gene from the same locus
- tag: C-terminal Strep-II
- host: E. coli BL21(DE3)
- purpose: co-expression and pull-down; test whether it is required for activity

### pLOC-full
- insert: the whole locus including the 21-copy array (23 bp repeats), native promoter
- tag: none
- host: E. coli MG1655 and a permissive native-adjacent host
- purpose: test array transcription and processing in context

### pLOC-array-deleted
- insert: the same locus with the array removed
- tag: none
- host: as above
- purpose: the array's contribution is only interpretable against its absence

## Assays

### RNA-dependent DNA polymerase assay
- readout: primer extension on a labelled RNA template, denaturing PAGE
- conditions: Mg2+ and Mn2+ titration, 30 and 37 C, +/- RNase H
- decides: whether the protein is an active RT at all

### Size-exclusion chromatography with multi-angle light scattering
- readout: oligomeric state of the enzyme alone and with the partner
- conditions: physiological salt
- decides: whether enzyme and partner form a defined complex

### Small-RNA sequencing and Northern blot
- readout: size and boundaries of transcripts from the array
- conditions: native host if available, otherwise the cloned locus
- decides: whether the array is expressed as discrete short RNAs, as the proposal requires

### Co-purifying nucleic acid sequencing
- readout: identity, strand and endpoints of nucleic acid bound to the enzyme
- conditions: pull-down from the host carrying the full locus
- decides: whether the array is the template the enzyme copies

### Efficiency of plaquing
- readout: plaque counts on hosts carrying the locus versus empty vector
- conditions: panel of coliphages, 30 and 37 C
- decides: whether the locus has a phage-defence phenotype

### Structure determination
- readout: cryo-EM or crystallography of the enzyme, and of the complex
- conditions: with and without the partner; with nucleic acid where it co-purifies
- decides: how the partner and the nucleic acid are positioned relative to the active site

## Controls

- empty vector in the same host and induction conditions
- catalytically dead active-site mutant for every activity readout
- a described RT from the catalogue (a retron RT) as a positive control for the polymerase assay
- mock pull-down from a host not expressing the tagged protein

## Scope and biosafety

- Standard BSL-1 work: heterologous expression of a protein of unknown function in a laboratory E. coli strain, and characterisation in vitro.
- Phage-challenge experiments use established laboratory coliphages and laboratory host strains.
- The handoff covers characterisation only. Any result suggesting a phenotype beyond the scope above should be reviewed before follow-up.

## Other loci with this architecture

Independent instances supporting the same proposal. They are evidence for the system, not additional experiments.

- `cand-ctg00187_g001`
- `cand-ctg00184_g001`
- `cand-ctg00192_g001`
- `cand-ctg00180_g001`
- `cand-ctg00193_g001`
- `cand-ctg00190_g001`
- `cand-ctg00183_g002`
- `cand-ctg00182_g001`
- `cand-ctg00188_g001`
- `cand-ctg00186_g002`
- `cand-ctg00189_g001`
- `cand-ctg00181_g001`
- `cand-ctg00179_g001`
