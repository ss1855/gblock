# amoA_AOA — status: blocked, not run against DeQueen

## Why this gene target is not screened here

`gene_presence_screen.py`'s `GENE_ALIASES` for `amoA_AOB` and `amoA_AOA` are
identical (`["amoA", "ammonia monooxygenase subunit A"]`), and `build_query()`
does not add any taxonomic domain restriction (Bacteria vs. Archaea) to the
NCBI query — the only domain restriction in the whole pipeline comes from
which species list is passed in via `--tax-table`.

The DeQueen dataset (`tax_table_df.csv`) is 1583/1584 rows Domain=="Bacteria"
(the 1584th is an "Unclassified"/"Unknown" placeholder row, not Archaea).
Running `--gene amoA_AOA --tax-table tax_table_df.csv` would therefore not
test archaeal presence at all — it would just re-run the exact same
"amoA"/"ammonia monooxygenase subunit A" title search against the same
all-Bacteria organism list that `amoA_AOB` already covers, and reproduce
amoA_AOB's hit counts under a misleading `gene_target=amoA_AOA` label. That's
a category error, not a meaningful "zero hits because no Archaea are present"
result, so it was not run.

## What's actually needed

amoA_AOA requires a separate archaeal reference species list (not yet
sourced — see README) before `gene_presence_screen.py` can produce a
meaningful result for this target. Once that list exists, the existing
script works unmodified: `--gene amoA_AOA --tax-table <archaeal_tax_table.csv>`.

## Status: blocked pending archaeal reference source. Not a pipeline failure.
