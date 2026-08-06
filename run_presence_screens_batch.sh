#!/bin/bash
# Sequentially run gene_presence_screen.py for the 5 non-ureC, non-amoA_AOA
# gene targets. Sequential (not parallel) deliberately: NCBI rate limits are
# enforced per source IP server-side, so N concurrent processes each doing
# their own ~3req/s internal sleep would sum to N*3req/s against the same
# IP and risk 429s/soft-blocking across all of them at once. Each individual
# gene_presence_screen.py call is already resumable (skips species already
# in its --out file), so this wrapper is safe to kill and re-run.

set -uo pipefail
cd /home/user/gblock

EMAIL="sujan.sth1991@gmail.com"
GENES=(amoA_AOB nirK nirS nosZ_cladeI nosZ_cladeII)

for gene in "${GENES[@]}"; do
    out="${gene}_presence.csv"
    log="${gene}_presence_run.log"
    echo "=== [$(date -u +%FT%TZ)] Starting/resuming $gene -> $out ===" | tee -a "$log"
    python3 gene_presence_screen.py --gene "$gene" --tax-table tax_table_df.csv \
        --email "$EMAIL" --out "$out" >> "$log" 2>&1
    echo "=== [$(date -u +%FT%TZ)] Finished $gene (exit $?) ===" | tee -a "$log"
done

echo "=== [$(date -u +%FT%TZ)] ALL 5 GENE PRESENCE SCREENS COMPLETE ===" | tee -a batch_presence_screens.log
