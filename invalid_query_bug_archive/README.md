# Invalid nirK/nirS presence screen runs -- archived, not used

These files were produced by the ORIGINAL gene_presence_screen.py, before
a confirmed NCBI esearch backend bug was found and fixed (see git history /
ANCHOR_TERMS docstring in gene_presence_screen.py for the full writeup).

Bug: a [Title] phrase of 3+ words ending in "nitrite reductase" combined
with `AND "<organism>"[Organism]` reliably returned Count=0, even against
verified exact-title matches. This silently broke nirK's and nirS's
compound-phrase aliases across all 1583 species -- nirK_presence.csv here
shows 0/1583 gene-positive, and the partial nirS run (131/1583 species
before it was caught and killed) shows 0/131.

Both genes were re-screened from scratch with the client-side-verification
fix (search_gene_verified() in gene_presence_screen.py). Do not use these
archived files for anything downstream.
