#!/bin/bash
# Sequentially run sequence_retrieval.py for the 5 non-ureC, non-amoA_AOA
# gene targets' gene-positive species. Sequential for the same NCBI
# per-IP-rate-limit reason as run_presence_screens_batch.sh. Each individual
# sequence_retrieval.py call is resumable.

set -uo pipefail
cd /home/user/gblock

EMAIL="sujan.sth1991@gmail.com"
GENES=(amoA_AOB nirK nirS nosZ_cladeI nosZ_cladeII)

for gene in "${GENES[@]}"; do
    fasta="${gene}_sequences.fasta"
    log_csv="${gene}_retrieval_log.csv"
    run_log="${gene}_retrieval_run.log"
    echo "=== [$(date -u +%FT%TZ)] Starting/resuming $gene sequence retrieval ===" | tee -a "$run_log"
    python3 sequence_retrieval.py --presence-csv "${gene}_presence.csv" \
        --email "$EMAIL" --out-fasta "$fasta" --out-log "$log_csv" >> "$run_log" 2>&1
    echo "=== [$(date -u +%FT%TZ)] Finished $gene sequence retrieval (exit $?) ===" | tee -a "$run_log"
done

echo "=== [$(date -u +%FT%TZ)] ALL 5 SEQUENCE RETRIEVALS COMPLETE ===" | tee -a batch_sequence_retrieval.log
