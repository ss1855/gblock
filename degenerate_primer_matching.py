#!/usr/bin/env python3
"""
In silico primer matching for DEGENERATE primers (contain IUPAC ambiguity
codes: R,Y,S,W,K,M,B,D,H,V,N,I). This is the degenerate-primer counterpart
to primer_matching.py, which is specifically documented as non-degenerate
(ureC only). All 6 remaining gene targets (amoA_AOB, nirK, nirS,
nosZ_cladeI, nosZ_cladeII, and amoA_AOA once unblocked) use degenerate
primers, so mismatch counting must treat an IUPAC code at a primer position
as a match against ANY of its represented literal bases, not as a literal
character comparison.

IUPAC HANDLING
    A primer position with code X matches a real sequence base Y (always
    literal A/C/G/T in retrieved sequences) if Y is in IUPAC_BASES[X].
    'I' (inosine) is treated as matching any base, same as N, since it is
    conventionally used in primer design exactly for that purpose. Reverse
    complementing a degenerate primer complements the base-identity of each
    ambiguity code (e.g. R=A/G complements to Y=C/T), not a literal-string
    complement -- see IUPAC_COMPLEMENT.

Otherwise structurally identical to primer_matching.py: sliding-window
best match per primer, boundary-irregular/fragment detection, and a
length-plausibility filter to reject off-target matches once mismatch
tolerance is loosened (see primer_matching.py's LENGTH PLAUSIBILITY CHECK
docstring for why this matters).

USAGE
    python degenerate_primer_matching.py --fasta nirK_sequences.fasta \
        --fwd-primer ATYGGCGGVAYGGCGA --rev-primer GCCTCGATCAGRTTRTGGT \
        --max-mismatches 2 \
        --out-fasta nirK_amplicons.fasta \
        --out-report nirK_primer_match_report.csv
"""

import argparse
import csv

IUPAC_BASES = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"),
    "B": set("CGT"), "D": set("AGT"), "H": set("ACT"), "V": set("ACG"),
    "N": set("ACGT"), "I": set("ACGT"),  # inosine: treated as matching any base
}

IUPAC_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N", "I": "I",
}


def revcomp_primer(primer: str) -> str:
    return "".join(IUPAC_COMPLEMENT[b] for b in reversed(primer.upper()))


def base_matches(seq_base: str, primer_code: str) -> bool:
    return seq_base in IUPAC_BASES.get(primer_code, set(seq_base))


def best_match(seq: str, primer: str):
    """Slide `primer` (may contain IUPAC codes) across `seq` (literal
    ACGT/N), return (start, mismatches) for the lowest-mismatch position."""
    plen = len(primer)
    if plen > len(seq):
        return None, None
    best_pos, best_mm = None, plen + 1
    for i in range(len(seq) - plen + 1):
        window = seq[i:i + plen]
        mm = sum(1 for a, b in zip(window, primer) if not base_matches(a, b))
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
            seqs[header].append(line.upper())
    return {h: "".join(chunks) for h, chunks in seqs.items()}


def boundary_irregular_species(seqs, short_fraction=0.9):
    lengths = sorted(len(s) for s in seqs.values())
    median_len = lengths[len(lengths) // 2] if lengths else 0
    flagged = set()
    for header, seq in seqs.items():
        frame_ok = len(seq) % 3 == 0
        start_ok = seq[:3] in ("ATG", "GTG", "TTG")
        stop_ok = seq[-3:] in ("TAA", "TAG", "TGA")
        if not (frame_ok and start_ok and stop_ok):
            flagged.add(header)
    short_fragments = {h for h in flagged if len(seqs[h]) < short_fraction * median_len}
    return flagged, short_fragments, median_len


def main():
    parser = argparse.ArgumentParser(
        description="Locate degenerate (IUPAC) primer binding sites and "
                     "extract amplicons from retrieved nucleotide sequences"
    )
    parser.add_argument("--fasta", required=True, help="sequence_retrieval.py FASTA output")
    parser.add_argument("--fwd-primer", required=True)
    parser.add_argument("--rev-primer", required=True)
    parser.add_argument("--max-mismatches", type=int, default=2)
    parser.add_argument("--length-tolerance", type=int, default=15,
                         help="Reject a match combination whose amplicon length "
                              "deviates from the strict-match (<=2 mismatch) "
                              "true length by more than this many bp")
    parser.add_argument("--short-fraction", type=float, default=0.9)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-report", required=True)
    args = parser.parse_args()

    seqs = parse_fasta(args.fasta)
    fwd_primer = args.fwd_primer.upper().replace(" ", "")
    rev_primer_rc = revcomp_primer(args.rev_primer.upper().replace(" ", ""))

    flagged, short_fragments, median_len = boundary_irregular_species(
        seqs, args.short_fraction)
    print(f"Total sequences: {len(seqs)}")
    print(f"Median sequence length: {median_len}")
    print(f"Boundary-irregular (frame/start/stop check): {len(flagged)}")
    print(f"  ...of which short fragments (< {args.short_fraction:.0%} of median): "
          f"{len(short_fragments)}")

    # Pass 1: compute both primers' best-match position/mismatch for every
    # species, independent of --max-mismatches.
    matches = {}
    for header, seq in seqs.items():
        fwd_pos, fwd_mm = best_match(seq, fwd_primer)
        rev_pos, rev_mm = best_match(seq, rev_primer_rc)
        matches[header] = (fwd_pos, fwd_mm, rev_pos, rev_mm)

    strict_lengths = []
    for header, (fwd_pos, fwd_mm, rev_pos, rev_mm) in matches.items():
        if (fwd_mm is not None and rev_mm is not None
                and fwd_mm <= 2 and rev_mm <= 2 and rev_pos > fwd_pos):
            strict_lengths.append(rev_pos + len(rev_primer_rc) - fwd_pos)
    true_length = max(set(strict_lengths), key=strict_lengths.count) if strict_lengths else None
    if true_length is not None:
        print(f"True amplicon length (from {len(strict_lengths)} strict matches): "
              f"{true_length}bp")
    else:
        print("No strict (<=2 mismatch) matches found -- length-plausibility "
              "filter disabled, all length-passing combinations accepted.")

    report_rows = []
    amplicons = []

    for header, seq in seqs.items():
        species, protein_acc, nuc_acc = header.split("|")
        fwd_pos, fwd_mm, rev_pos, rev_mm = matches[header]

        fwd_ok = fwd_mm is not None and fwd_mm <= args.max_mismatches
        rev_ok = rev_mm is not None and rev_mm <= args.max_mismatches
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
                status = "off_target_position"
            elif not fwd_ok and not rev_ok:
                status = "no_match_both"
            elif not fwd_ok:
                status = "no_match_fwd"
            elif not rev_ok:
                status = "no_match_rev"
            else:
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
