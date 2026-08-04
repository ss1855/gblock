#!/usr/bin/env python3
"""
Diagnostic (not a pipeline stage): characterize WHY Firmicutes ureC sequences
are failing primer matching, ahead of deciding whether the fix is a
Firmicutes-specific primer redesign, accepting the primer's real taxonomic
boundary, or a gene-presence-screen annotation gap (fused ureABC operons).

Unlike primer_matching.py, this reports the actual best-alignment mismatch
count/position for both primers regardless of threshold, plus a per-position
mismatch pattern oriented to each primer's own 5'->3' frame (so 3'-terminal
mismatches, which matter most for real PCR extension, are identifiable).

USAGE
    python diagnose_firmicutes_failures.py --fasta ureC_sequences.fasta \
        --presence-csv ureC_presence.csv --report-csv ureC_primer_match_report.csv \
        --fwd-primer TTCACACCTTCCACACCGAA --rev-primer AACGTCGGGTTGGTCGAG \
        --sample-size 18
"""

import argparse
import csv
import random


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


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


def mismatch_at_position(seq, primer, pos):
    """Mismatch count for `primer` against the fixed window seq[pos:pos+len(primer)]."""
    if pos is None or pos + len(primer) > len(seq):
        return None
    window = seq[pos:pos + len(primer)]
    return sum(1 for a, b in zip(window, primer) if a != b)


def primer_oriented_mismatches(primer, seq, pos, is_reverse):
    """Return a list of booleans, index 0 = primer's 5' end, aligned to
    `primer`'s own 5'->3' orientation regardless of fwd/rev."""
    window = seq[pos:pos + len(primer)]
    aligned = revcomp(window) if is_reverse else window
    return [a != b for a, b in zip(primer, aligned)]


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--presence-csv", required=True)
    parser.add_argument("--report-csv", required=True)
    parser.add_argument("--fwd-primer", required=True)
    parser.add_argument("--rev-primer", required=True)
    parser.add_argument("--sample-size", type=int, default=18)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    seqs = parse_fasta(args.fasta)
    fwd_primer = args.fwd_primer.upper()
    rev_primer = args.rev_primer.upper()

    presence = {}
    with open(args.presence_csv, newline="") as f:
        for row in csv.DictReader(f):
            presence[row["species"]] = row

    status_by_species = {}
    with open(args.report_csv, newline="") as f:
        for row in csv.DictReader(f):
            status_by_species[row["species"]] = row["status"]

    # Establish the "expected" homologous position from the primer matching
    # report's own true-length reference isn't available here, so re-derive
    # it the same way: strict (<=2 mismatch) matches across ALL species.
    strict_fwd_pos, strict_rev_pos = [], []
    for header, seq in seqs.items():
        fp, fm = best_match(seq, fwd_primer)
        rp, rm = best_match(seq, revcomp(rev_primer))
        if fm is not None and fm <= 2:
            strict_fwd_pos.append(fp)
        if rm is not None and rm <= 2:
            strict_rev_pos.append(rp)
    expected_fwd_pos = sorted(strict_fwd_pos)[len(strict_fwd_pos) // 2]
    expected_rev_pos = sorted(strict_rev_pos)[len(strict_rev_pos) // 2]
    print(f"Expected fwd primer position (median of strict matches): {expected_fwd_pos}")
    print(f"Expected rev primer position (median of strict matches): {expected_rev_pos}")
    print()

    # Firmicutes species present in the FASTA that did not succeed.
    failing = []
    for header, seq in seqs.items():
        species = header.split("|")[0]
        row = presence.get(species)
        if not row or row["phylum"] != "Firmicutes":
            continue
        status = status_by_species.get(species, "ABSENT_FROM_REPORT")
        if status != "success":
            failing.append((species, row["genus"], header, seq, status))

    print(f"Failing Firmicutes present in FASTA: {len(failing)}")
    genera = sorted(set(g for _, g, _, _, _ in failing))
    print(f"Distinct genera: {len(genera)}")

    random.seed(args.seed)
    random.shuffle(genera)
    chosen_genera = genera[:args.sample_size]
    sample = []
    for g in chosen_genera:
        candidates = [f for f in failing if f[1] == g]
        sample.append(random.choice(candidates))
    sample.sort(key=lambda x: x[1])

    print(f"\nSampled {len(sample)} species across {len(chosen_genera)} genera:\n")

    results = []
    for species, genus, header, seq, status in sample:
        fwd_pos, fwd_mm = best_match(seq, fwd_primer)
        rev_pos, rev_mm = best_match(seq, revcomp(rev_primer))

        fwd_at_expected = mismatch_at_position(seq, fwd_primer, expected_fwd_pos)
        rev_at_expected = mismatch_at_position(seq, revcomp(rev_primer), expected_rev_pos)

        fwd_mm_pattern = primer_oriented_mismatches(fwd_primer, seq, fwd_pos, is_reverse=False)
        rev_mm_pattern = primer_oriented_mismatches(rev_primer, seq, rev_pos, is_reverse=True)

        def three_prime_mm(pattern, n=5):
            return sum(pattern[-n:])

        # Categorize: is the globally-best-scoring window both low-mismatch
        # (well below ~75% random expectation for 4-letter DNA) AND within a
        # generous window of the known homologous position (small indels
        # elsewhere in the sequence can shift a pure-Hamming local window by
        # a few bp without indicating the site is a different locus)? If
        # so, the site is present but diverged. If the sequence is too
        # short to even reach the expected position, that's a fragment
        # issue, not absence. Otherwise, treat as genuinely absent.
        POS_WINDOW = 15
        MM_CEILING = 8  # well below ~15-16/20 expected from a random window
        fwd_near_expected = fwd_pos is not None and abs(fwd_pos - expected_fwd_pos) <= POS_WINDOW
        rev_near_expected = rev_pos is not None and abs(rev_pos - expected_rev_pos) <= POS_WINDOW
        fwd_plausible = fwd_mm is not None and fwd_mm <= MM_CEILING
        rev_plausible = rev_mm is not None and rev_mm <= MM_CEILING

        fwd_site_reachable = expected_fwd_pos + len(fwd_primer) <= len(seq)
        rev_site_reachable = expected_rev_pos + len(rev_primer) <= len(seq)

        if not fwd_site_reachable and not rev_site_reachable:
            category = "fragment_site_unreachable"
        elif (fwd_near_expected and fwd_plausible) or (rev_near_expected and rev_plausible):
            category = "diverged_present"
        else:
            category = "genuinely_absent"

        results.append({
            "species": species, "genus": genus, "status": status,
            "seq_len": len(seq),
            "fwd_best_pos": fwd_pos, "fwd_best_mm": fwd_mm,
            "fwd_mm_at_expected_pos": fwd_at_expected,
            "fwd_near_expected_pos": fwd_near_expected,
            "fwd_3prime_mm_of5": three_prime_mm(fwd_mm_pattern) if fwd_pos is not None else None,
            "rev_best_pos": rev_pos, "rev_best_mm": rev_mm,
            "rev_mm_at_expected_pos": rev_at_expected,
            "rev_near_expected_pos": rev_near_expected,
            "rev_3prime_mm_of5": three_prime_mm(rev_mm_pattern) if rev_pos is not None else None,
            "category": category,
        })

        print(f"{species} ({genus}) [{status}] seq_len={len(seq)}")
        print(f"  FWD best: pos={fwd_pos} mm={fwd_mm}/{len(fwd_primer)}  "
              f"(at expected pos {expected_fwd_pos}: mm={fwd_at_expected})  "
              f"3'-end(5bp) mismatches={three_prime_mm(fwd_mm_pattern) if fwd_pos is not None else 'NA'}")
        print(f"  REV best: pos={rev_pos} mm={rev_mm}/{len(rev_primer)}  "
              f"(at expected pos {expected_rev_pos}: mm={rev_at_expected})  "
              f"3'-end(5bp) mismatches={three_prime_mm(rev_mm_pattern) if rev_pos is not None else 'NA'}")
        print(f"  => category: {category}")
        print()

    with open("firmicutes_diagnostic_sample.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    print("=== Category summary (sample of", len(results), ") ===")
    from collections import Counter
    print(Counter(r["category"] for r in results))
    print("\nWritten: firmicutes_diagnostic_sample.csv")

    # Full-scale pass (not just the sample) for statistically reliable
    # category fractions across every failing Firmicutes species.
    print(f"\n=== Full-scale categorization across all {len(failing)} failing "
          f"Firmicutes ===")
    full_results = []
    for species, genus, header, seq, status in failing:
        fwd_pos, fwd_mm = best_match(seq, fwd_primer)
        rev_pos, rev_mm = best_match(seq, revcomp(rev_primer))

        fwd_near_expected = fwd_pos is not None and abs(fwd_pos - expected_fwd_pos) <= POS_WINDOW
        rev_near_expected = rev_pos is not None and abs(rev_pos - expected_rev_pos) <= POS_WINDOW
        fwd_plausible = fwd_mm is not None and fwd_mm <= MM_CEILING
        rev_plausible = rev_mm is not None and rev_mm <= MM_CEILING
        fwd_site_reachable = expected_fwd_pos + len(fwd_primer) <= len(seq)
        rev_site_reachable = expected_rev_pos + len(rev_primer) <= len(seq)

        if not fwd_site_reachable and not rev_site_reachable:
            category = "fragment_site_unreachable"
        elif (fwd_near_expected and fwd_plausible) or (rev_near_expected and rev_plausible):
            category = "diverged_present"
        else:
            category = "genuinely_absent"

        full_results.append({
            "species": species, "genus": genus, "status": status,
            "seq_len": len(seq), "fwd_best_pos": fwd_pos, "fwd_best_mm": fwd_mm,
            "rev_best_pos": rev_pos, "rev_best_mm": rev_mm, "category": category,
        })

    with open("firmicutes_diagnostic_full.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(full_results[0].keys()))
        writer.writeheader()
        writer.writerows(full_results)

    counts = Counter(r["category"] for r in full_results)
    total = len(full_results)
    for cat, n in counts.most_common():
        print(f"  {cat}: {n} ({100*n/total:.1f}%)")
    print("\nWritten: firmicutes_diagnostic_full.csv")

    # Per-position mismatch frequency across all diverged_present species,
    # to check whether mismatches concentrate at primer 3' ends (real PCR
    # failure signature) or follow a codon wobble-position pattern
    # (synonymous substitutions -- protein-conserved, DNA-degenerate).
    print("\n=== Per-position mismatch frequency (diverged_present only) ===")
    diverged_headers = [(species, header, seq) for species, genus, header, seq, status in failing
                         for r in full_results
                         if r["species"] == species and r["category"] == "diverged_present"]
    fwd_pos_mm = [0] * len(fwd_primer)
    rev_pos_mm = [0] * len(rev_primer)
    n_div = 0
    seen = set()
    for species, header, seq in diverged_headers:
        if species in seen:
            continue
        seen.add(species)
        n_div += 1
        fp, fm = best_match(seq, fwd_primer)
        rp, rm = best_match(seq, revcomp(rev_primer))
        window = seq[fp:fp + len(fwd_primer)]
        for i, (a, b) in enumerate(zip(fwd_primer, window)):
            if a != b:
                fwd_pos_mm[i] += 1
        rwindow = revcomp(seq[rp:rp + len(rev_primer)])
        for i, (a, b) in enumerate(zip(rev_primer, rwindow)):
            if a != b:
                rev_pos_mm[i] += 1

    print(f"n={n_div}")
    print(f"FWD primer 5'->3': {fwd_primer}")
    print("  position: mismatch%% -> " +
          ", ".join(f"{i+1}:{100*c/n_div:.0f}%%" for i, c in enumerate(fwd_pos_mm)))
    print(f"REV primer 5'->3': {rev_primer}")
    print("  position: mismatch%% -> " +
          ", ".join(f"{i+1}:{100*c/n_div:.0f}%%" for i, c in enumerate(rev_pos_mm)))


if __name__ == "__main__":
    main()
