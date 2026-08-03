# gblock
# gBlock Consensus Amplicon Design (Pilot)

## Project overview

Designing gBlock positive controls for a digital PCR (dPCR) assay quantifying
nitrogen cycle genes in litter-dwelling bacteria. Because the primers are
degenerate, a single gBlock per gene needs to be built from a **consensus
amplicon sequence** representative of the target bacterial community, rather
than an arbitrary single-organism sequence.

The [DeQueen dataset](./tax_table_df.csv) (1584 bacterial species, identified
to species level, spanning poultry litter-relevant phyla) is being used to
guide this. Full taxonomy: Firmicutes (808), Actinobacteria (503),
Proteobacteria (202), Bacteroidota (54), plus several minor phyla.

## Gene targets

Degenerate primers for 7 targets, spanning the broader nitrogen cycle
(ammonification, nitrification, denitrification), not just ammonia
biosynthesis narrowly. See
[`IDT_degenerate_primer_063026.xlsx`](./IDT_degenerate_primer_063026.xlsx)
for the full primer sequences.

| Gene | Process | Notes |
|---|---|---|
| ureC | Ammonification (urea hydrolysis) | Broadly distributed |
| amoA (AOB) | Nitrification (bacterial) | Narrow phylogenetic distribution (Proteobacteria) |
| amoA (AOA) | Nitrification (archaeal) | **Targets Archaea, not Bacteria** — DeQueen (all Bacteria) can't supply references for this one; needs a separate archaeal reference set |
| nirK | Denitrification | Narrow distribution |
| nirS | Denitrification | Narrow distribution |
| nosZ clade I | Denitrification | Narrow distribution |
| nosZ clade II | Denitrification | Narrow distribution, uses inosine in primer (can't be directly synthesized into a gBlock template) |

## Key methodological point

These are **functional marker genes, not universal genes like 16S**. Most of
the 1584 DeQueen species will not carry any single given target. The
pipeline therefore screens for gene presence first, then builds the
consensus only from the gene-positive subset, rather than assuming
every species carries every gene.

## Pipeline stages

1. **Gene presence screen** (`gene_presence_screen.py`) — queries NCBI's
   protein database per species per gene to determine which DeQueen taxa
   actually carry each target. Resumable, gene selectable via `--gene` flag.
2. **Sequence retrieval** — pull actual gene sequences for the gene-positive
   subset per target (not yet implemented).
3. **In silico primer matching** — expand degenerate primer positions,
   extract amplicon region from each gene-positive sequence, report
   per-taxon match rates (not yet implemented).
4. **Alignment and consensus** — MAFFT/MUSCLE align extracted amplicons,
   build majority-rule or medoid consensus per gene (not yet implemented).
5. **Validation** — back-check the finished consensus against the primers
   to confirm no mismatch was introduced at the flanks (not yet implemented).

## Status

Pilot phase: validating the gene presence screen end to end on `ureC`
before running the full panel across all 7 targets.

## Usage

```bash
pip install biopython pandas

# Dry run (no network calls, sanity checks query construction)
python gene_presence_screen.py --gene ureC --tax-table tax_table_df.csv \
    --email your_email@unt.edu --dry-run --limit 5

# Small live pilot
python gene_presence_screen.py --gene ureC --tax-table tax_table_df.csv \
    --email your_email@unt.edu --api-key YOUR_NCBI_API_KEY \
    --out ureC_pilot.csv --limit 20

# Full run (resumable)
python gene_presence_screen.py --gene ureC --tax-table tax_table_df.csv \
    --email your_email@unt.edu --api-key YOUR_NCBI_API_KEY \
    --out ureC_presence.csv
```

## Requirements

- Outbound internet access to `eutils.ncbi.nlm.nih.gov` (required for the
  gene presence screen — if running via Claude Code on the web, set the
  environment's network access to full internet access, the default
  allowlisted mode does not include NCBI)
- NCBI account API key recommended (raises rate limit from ~3 to ~10 req/s)

## Collaborators

Sujan (bioinformatics), Gaurav, Brett Hale (AgriGro, Head of R&D)
