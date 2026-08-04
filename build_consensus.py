#!/usr/bin/env python3
"""
Alignment and consensus (pipeline stage 4): from a MAFFT-aligned amplicon
FASTA, compute per-column conservation/entropy, a gap-aware majority-rule
consensus, and the medoid sequence (the real species sequence with minimum
average distance to all others). Reports how much the two candidate
representative sequences differ, and back-checks both against the original
primers to confirm the consensus-building step didn't introduce a mismatch
at the primer flanks.

USAGE
    python build_consensus.py --aligned-fasta ureC_amplicons_aligned.fasta \
        --fwd-primer TTCACACCTTCCACACCGAA --rev-primer AACGTCGGGTTGGTCGAG \
        --max-mismatches 4 \
        --out-consensus ureC_consensus.fasta \
        --out-conservation ureC_conservation.csv

REQUIREMENTS
    none beyond the standard library
"""

import argparse
import csv
import math
from collections import Counter


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTN-", "TGCAN-")
    return seq.translate(comp)[::-1]


def parse_fasta(path):
    seqs = {}
    header = None
    for line in open(path):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:]
            seqs[header] = []
        else:
            seqs[header].append(line)
    return {h: "".join(c) for h, c in seqs.items()}


def column_entropy(column):
    counts = Counter(column)
    n = len(column)
    ent = 0.0
    for c, n_c in counts.items():
        p = n_c / n
        ent -= p * math.log2(p)
    return ent


def majority_consensus(seqs_list, aln_len):
    """Gap-aware majority vote per column. Ties broken by first-seen base
    in a fixed priority order (A,C,G,T,-) for determinism."""
    consensus = []
    priority = "ACGT-"
    for col in range(aln_len):
        counts = Counter(s[col] for s in seqs_list)
        best_base = max(counts.items(), key=lambda kv: (kv[1], -priority.index(kv[0])))[0]
        consensus.append(best_base)
    return "".join(consensus).replace("-", "")


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def find_medoid(headers, seqs_list):
    """Real sequence with minimum total Hamming distance to all others."""
    n = len(seqs_list)
    totals = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            d = hamming(seqs_list[i], seqs_list[j])
            totals[i] += d
            totals[j] += d
    best_i = min(range(n), key=lambda i: totals[i])
    return headers[best_i], seqs_list[best_i], totals[best_i] / (n - 1)


def best_match(seq, primer):
    plen = len(primer)
    if plen > len(seq):
        return None, None
    best_pos, best_mm = None, plen + 1
    for i in range(len(seq) - plen + 1):
        window = seq[i:i + plen]
        mm = sum(1 for a, b in zip(window, primer) if a != b)
        if mm < best_mm:
            best_pos, best_mm = i, mm
    return best_pos, best_mm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-fasta", required=True)
    parser.add_argument("--fwd-primer", required=True)
    parser.add_argument("--rev-primer", required=True)
    parser.add_argument("--max-mismatches", type=int, default=4,
                         help="Threshold to report pass/fail on the primer back-check")
    parser.add_argument("--out-consensus", required=True)
    parser.add_argument("--out-conservation", required=True)
    args = parser.parse_args()

    seqs = parse_fasta(args.aligned_fasta)
    headers = list(seqs.keys())
    seqs_list = [seqs[h].upper() for h in headers]
    aln_len = len(seqs_list[0])
    n = len(seqs_list)
    print(f"Alignment: {n} sequences x {aln_len} columns")

    # --- Per-column conservation / entropy ---
    conservation_rows = []
    split_columns = []
    for col in range(aln_len):
        column = [s[col] for s in seqs_list]
        counts = Counter(column)
        ent = column_entropy(column)
        top_base, top_n = counts.most_common(1)[0]
        top_frac = top_n / n
        # "genuinely split" -- not just one dominant base with rare noise:
        # second-most-common base present in >10% of sequences.
        second_frac = (sorted(counts.values(), reverse=True)[1] / n) if len(counts) > 1 else 0.0
        is_split = second_frac > 0.10
        conservation_rows.append({
            "column": col + 1, "entropy_bits": round(ent, 4),
            "top_base": top_base, "top_base_frac": round(top_frac, 4),
            "second_base_frac": round(second_frac, 4),
            "distinct_bases": len(counts), "genuinely_split": is_split,
        })
        if is_split:
            split_columns.append(col + 1)

    with open(args.out_conservation, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(conservation_rows[0].keys()))
        writer.writeheader()
        writer.writerows(conservation_rows)

    print(f"Genuinely split columns (second allele >10% frequency): "
          f"{len(split_columns)}/{aln_len} -> {split_columns}")

    # --- Majority-rule consensus (gap-aware) ---
    majority = majority_consensus(seqs_list, aln_len)
    print(f"\nMajority-rule consensus ({len(majority)}bp): {majority}")

    # --- Medoid ---
    medoid_header, medoid_seq_aligned, medoid_avg_dist = find_medoid(headers, seqs_list)
    medoid_seq = medoid_seq_aligned.replace("-", "")
    print(f"\nMedoid: {medoid_header}")
    print(f"  avg Hamming distance to all others: {medoid_avg_dist:.3f}")
    print(f"  sequence ({len(medoid_seq)}bp): {medoid_seq}")

    # --- Compare majority-rule vs medoid ---
    if len(majority) == len(medoid_seq):
        diff = hamming(majority, medoid_seq)
        pct_identity = 100 * (1 - diff / len(majority))
        print(f"\nMajority-rule vs medoid: {diff} differences / {len(majority)}bp "
              f"({pct_identity:.1f}% identity)")
    else:
        diff, pct_identity = None, None
        print(f"\nMajority-rule ({len(majority)}bp) and medoid ({len(medoid_seq)}bp) "
              f"differ in length -- compare via alignment, not raw Hamming.")

    # --- Back-check both candidates against the original primers ---
    fwd_primer = args.fwd_primer.upper()
    rev_primer_rc = revcomp(args.rev_primer.upper())
    print("\n=== Primer back-check ===")
    for name, seq in [("majority-rule consensus", majority), ("medoid", medoid_seq)]:
        fp, fmm = best_match(seq, fwd_primer)
        rp, rmm = best_match(seq, rev_primer_rc)
        fwd_ok = fmm is not None and fmm <= args.max_mismatches
        rev_ok = rmm is not None and rmm <= args.max_mismatches
        print(f"{name}: fwd pos={fp} mm={fmm}/{len(fwd_primer)} "
              f"({'OK' if fwd_ok else 'FAIL'}); "
              f"rev pos={rp} mm={rmm}/{len(args.rev_primer)} "
              f"({'OK' if rev_ok else 'FAIL'})")

    with open(args.out_consensus, "w") as f:
        f.write(f">ureC_majority_rule_consensus n={n}\n{majority}\n")
        f.write(f">ureC_medoid|{medoid_header}\n{medoid_seq}\n")

    print(f"\nWritten: {args.out_consensus}, {args.out_conservation}")


if __name__ == "__main__":
    main()
