#!/usr/bin/env python3
"""
In silico primer matching: locate the ureC (or any non-degenerate primer
pair's) binding sites in each retrieved nucleotide sequence and extract the
amplicon region between them.

WHY THIS STEP MATTERS
    sequence_retrieval.py pulls the full coding sequence per species, but the
    gBlock only needs the amplicon between the primer pair. This script finds
    each primer's best-matching position (allowing a small number of
    mismatches, since these are real species and not the primer's design
    organism), reports match quality per species, and extracts the amplicon
    for species that pass. It is local sequence analysis only -- no network
    calls, so no rate limiting/resumability machinery is needed.

MATCHING
    - Forward primer is searched directly against the sequence (sense
      strand). Reverse primer is searched as its reverse complement, since
      retrieved sequences are all in coding/sense-strand orientation
      (sequence_retrieval.py already reverse-complements minus-strand CDS
      regions). Both searches use a sliding-window Hamming distance (no
      indels) and keep the position with the fewest mismatches.
    - A primer is considered "matched" if its best position has
      <= --max-mismatches (default 2) mismatches. The best mismatch count is
      always reported, even above the threshold, since a 3-4 mismatch
      near-miss is more informative than a bare "no match."

LENGTH PLAUSIBILITY CHECK
    Loosening --max-mismatches to capture more species opens a failure mode
    the strict (<=2 mismatch) threshold doesn't have: the forward and
    reverse primer's independently-best-scoring windows can each be within
    the mismatch threshold while sitting at completely different, unrelated
    positions in the sequence (a coincidental N-mismatch match far from the
    true conserved locus, rather than the true homologous site which may
    have had one extra mismatch and lost the "best window" comparison). This
    produces amplicons with implausible lengths. To catch this, the true
    amplicon length is established from the strict (<=2 mismatch) matches
    only -- those are unambiguous -- and any match (at any mismatch
    threshold) whose amplicon length falls outside
    [true_length - --length-tolerance, true_length + --length-tolerance] is
    rejected as "off_target_position" rather than accepted.

BOUNDARY-IRREGULAR SEQUENCES
    Before matching, every sequence is checked for frame (length % 3 == 0)
    and canonical start/stop codons. Flagged sequences that are also
    substantially shorter than the typical full-length CDS
    (< --short-fraction of the median retrieved length) are treated as
    genuine fragments: if primer matching then fails for one of these, it is
    logged as "insufficient_sequence" rather than "no_match_*", since the
    failure is attributable to a truncated record, not primer/species
    divergence. Boundary-irregular sequences that are near full length
    (e.g. an unusual start codon on an otherwise complete CDS) are matched
    normally -- they are not assumed to be fragments.

USAGE
    python primer_matching.py --fasta ureC_sequences.fasta \
        --fwd-primer TTCACACCTTCCACACCGAA --rev-primer AACGTCGGGTTGGTCGAG \
        --out-fasta ureC_amplicons.fasta \
        --out-report ureC_primer_match_report.csv

REQUIREMENTS
    none beyond the standard library
"""

import argparse
import csv


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def best_match(seq: str, primer: str):
    """Slide `primer` across `seq`, return (start, mismatches) for the best
    (lowest-mismatch) position, or (None, None) if primer is longer than seq."""
    plen = len(primer)
    if plen > len(seq):
        return None, None
    best_pos, best_mm = None, plen + 1
    for i in range(len(seq) - plen + 1):
        window = seq[i:i + plen]
        mm = sum(1 for a, b in zip(window, primer) if a != b)
        if mm < best_mm:
            best_pos, best_mm = i, mm
            if mm == 0:
                break
    return best_pos, best_mm


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
    return {h: "".join(chunks) for h, chunks in seqs.items()}


def boundary_irregular_species(seqs, short_fraction=0.9):
    lengths = sorted(len(s) for s in seqs.values())
    median_len = lengths[len(lengths) // 2]
    flagged = set()
    for header, seq in seqs.items():
        frame_ok = len(seq) % 3 == 0
        start_ok = seq[:3] in ("ATG", "GTG", "TTG")
        stop_ok = seq[-3:] in ("TAA", "TAG", "TGA")
        if not (frame_ok and start_ok and stop_ok):
            flagged.add(header)
    short_fragments = {h for h in flagged if len(seqs[h]) < short_fraction * median_len}
    return flagged, short_fragments, median_len


def classify_mismatch(mm, max_mismatches):
    if mm is None:
        return "no_window"
    if mm == 0:
        return "perfect"
    if mm <= max_mismatches:
        return f"{mm}_mismatch"
    return "no_match"


def main():
    parser = argparse.ArgumentParser(
        description="Locate primer binding sites and extract amplicons from "
                     "retrieved nucleotide sequences"
    )
    parser.add_argument("--fasta", required=True, help="sequence_retrieval.py FASTA output")
    parser.add_argument("--fwd-primer", required=True)
    parser.add_argument("--rev-primer", required=True)
    parser.add_argument("--max-mismatches", type=int, default=2)
    parser.add_argument("--length-tolerance", type=int, default=15,
                         help="Reject a match combination whose amplicon length "
                              "deviates from the strict-match (<=2 mismatch) "
                              "true length by more than this many bp")
    parser.add_argument("--short-fraction", type=float, default=0.9,
                         help="Boundary-irregular sequences shorter than this "
                              "fraction of the median length are treated as "
                              "fragments for insufficient_sequence logging")
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-report", required=True)
    args = parser.parse_args()

    seqs = parse_fasta(args.fasta)
    fwd_primer = args.fwd_primer.upper()
    rev_primer_rc = revcomp(args.rev_primer.upper())

    flagged, short_fragments, median_len = boundary_irregular_species(
        seqs, args.short_fraction)
    print(f"Median sequence length: {median_len}")
    print(f"Boundary-irregular (frame/start/stop check): {len(flagged)}")
    print(f"  ...of which short fragments (< {args.short_fraction:.0%} of median): "
          f"{len(short_fragments)}")

    # Pass 1: compute both primers' best-match position/mismatch count for
    # every species, independent of --max-mismatches.
    matches = {}
    for header, seq in seqs.items():
        fwd_pos, fwd_mm = best_match(seq, fwd_primer)
        rev_pos, rev_mm = best_match(seq, rev_primer_rc)
        matches[header] = (fwd_pos, fwd_mm, rev_pos, rev_mm)

    # Establish the true amplicon length from strict (<=2 mismatch) matches
    # only -- those are unambiguous, so their length is a reliable reference
    # for rejecting off-target matches once the threshold is loosened.
    strict_lengths = []
    for header, (fwd_pos, fwd_mm, rev_pos, rev_mm) in matches.items():
        if (fwd_mm is not None and rev_mm is not None
                and fwd_mm <= 2 and rev_mm <= 2 and rev_pos > fwd_pos):
            strict_lengths.append(rev_pos + len(rev_primer_rc) - fwd_pos)
    true_length = max(set(strict_lengths), key=strict_lengths.count) if strict_lengths else None
    if true_length is not None:
        print(f"True amplicon length (from {len(strict_lengths)} strict matches): "
              f"{true_length}bp")

    report_rows = []
    amplicons = []

    for header, seq in seqs.items():
        species, protein_acc, nuc_acc = header.split("|")
        fwd_pos, fwd_mm, rev_pos, rev_mm = matches[header]

        fwd_ok = fwd_mm is not None and fwd_mm <= args.max_mismatches
        rev_ok = rev_mm is not None and rev_mm <= args.max_mismatches
        # Reverse primer's binding site must come after the forward primer's.
        order_ok = fwd_ok and rev_ok and rev_pos > fwd_pos

        amplicon_len = ""
        length_ok = True
        if order_ok:
            amplicon_len = rev_pos + len(rev_primer_rc) - fwd_pos
            if true_length is not None:
                length_ok = abs(amplicon_len - true_length) <= args.length_tolerance

        if order_ok and length_ok:
            amp_start = fwd_pos
            amp_end = rev_pos + len(rev_primer_rc)
            amplicon = seq[amp_start:amp_end]
            status = "success"
            amplicons.append((header, amplicon))
        else:
            if header in short_fragments:
                status = "insufficient_sequence"
            elif order_ok and not length_ok:
                # both primers matched within threshold, but at positions
                # implying an implausible amplicon length -- an off-target
                # coincidental match, not the true conserved locus.
                status = "off_target_position"
            elif not fwd_ok and not rev_ok:
                status = "no_match_both"
            elif not fwd_ok:
                status = "no_match_fwd"
            elif not rev_ok:
                status = "no_match_rev"
            else:
                # both matched individually but reverse primer sits before
                # forward primer -- inconsistent orientation, not usable.
                status = "no_match_order"

        report_rows.append({
            "species": species,
            "mismatch_count_fwd": fwd_mm if fwd_mm is not None else "",
            "mismatch_count_rev": rev_mm if rev_mm is not None else "",
            "amplicon_length": amplicon_len,
            "status": status,
            "boundary_irregular": header in flagged,
        })

    with open(args.out_report, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "species", "mismatch_count_fwd", "mismatch_count_rev",
            "amplicon_length", "status", "boundary_irregular"])
        writer.writeheader()
        for row in report_rows:
            writer.writerow(row)

    with open(args.out_fasta, "w") as f:
        for header, amplicon in amplicons:
            f.write(f">{header}\n{amplicon}\n")

    statuses = [r["status"] for r in report_rows]
    print(f"\nTotal species: {len(report_rows)}")
    for s in sorted(set(statuses)):
        print(f"  {s}: {statuses.count(s)}")
    print(f"\nAmplicons written to {args.out_fasta}")
    print(f"Report written to {args.out_report}")


if __name__ == "__main__":
    main()
