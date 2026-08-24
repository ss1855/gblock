# ureC "3b" regression diagnostic — impact analysis

## Source
`ureC_3b_regressions_diagnosed.csv` — an independent HMM+BLAST-based verification
(stage "3b", run outside this pipeline) cross-checking our original NCBI
title-text-based ureC presence screen. Flags 74 species where our screen called
ureC-positive but assembly-level HMM/BLAST search against a specific selected
genome assembly found no true ureC ortholog, with two diagnoses:
- "no urease-family gene at all in selected assembly" (different strain than
  our hit, not a false negative on 3b's part)
- "different urease-operon gene present (accessory/gamma subunit), not true
  ureC" (our hit matched a paralog, not true ureC)

## Validation
Cross-checked the file's "old_hit_count"/"old_best_accession" columns against
our actual `ureC_presence.csv` -- exact match, confirming this is a genuine
cross-reference against our real pipeline output, not unrelated/fabricated data.

Independently verified the E. coli case via direct NCBI query: our hit
(BHS50992.1) is from E. coli O157:H7; the standard reference strain (K-12
MG1655) has zero urease-titled protein records. This matches well-known
microbiology (E. coli is textbook urease-negative) and confirms the 3b
diagnostic's "different strain" explanation is accurate, not spurious.

## Scope
- 74 species flagged total; 64 (86%) are in our delivered 538-species
  consensus pool (12% of the pool).
- Roughly proportional across major phyla (Actinobacteria 8.6%, Firmicutes
  13.6%, Proteobacteria 12.9%); notably higher in Bacteroidota (55.6%, 5/9 --
  small subgroup).

## Impact on the delivered ureC consensus
Rebuilt the majority-rule consensus from the 538-pool alignment
(`ureC_msa_anchored_amplicons_aligned.fasta`) excluding all 64 flagged
species: **0/103 positions differ** from the delivered consensus
(`ureC_consensus_538.fasta`). The majority-rule method is robust to this
level of contamination -- 12% of the pool isn't enough weight to flip the
majority call at any position.

## Conclusion
The ureC gBlock candidate sequence does not need to be revised based on this
finding. However, it confirms a real, general methodological gap in the
title-text-based gene presence screen (species-level match does not verify
strain/assembly-level accuracy) that likely also affects the still-pending
gene targets (nirK, amoA_AOB, nosZ_cladeI, nosZ_cladeII). Worth reconciling
with whoever runs the "3a/3b" verification pipeline so both methods converge,
and worth considering for any future gene presence screening in this project.
