#!/usr/bin/env python3
"""
Protein-guided (codon-aware) nucleotide alignment. For gene families with
high nucleotide divergence but conserved protein sequence (synonymous
substitutions dominate), a raw nucleotide MAFFT alignment can fail to
correctly place homologous regions across a taxonomically diverse pool --
this was observed for nosZ_cladeI (a 267bp anchored locus produced a
1388-column, extremely gap-heavy alignment window, recovering only 1/184
species instead of the broad recovery seen for ureC).

This script translates each nucleotide sequence to protein, aligns the
PROTEINS with MAFFT (protein-level similarity is far more robust to
synonymous drift), then maps the protein alignment back to a codon-aligned
nucleotide alignment (each protein alignment column -> 3 nucleotide
columns, gaps -> "---"). The resulting nucleotide alignment can then be
used for the same anchor-mapping / MSA-anchored extraction as before.

USAGE
    python protein_guided_alignment.py --fasta nosZ_cladeI_sequences_full_pool.fasta \
        --out-protein-fasta nosZ_cladeI_proteins.fasta \
        --out-nt-alignment nosZ_cladeI_codon_aligned.fasta
"""

import argparse
import subprocess
import sys


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
    return {h: "".join(c) for h, c in seqs.items()}


CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def translate(seq):
    seq = seq[:len(seq) - len(seq) % 3]  # trim to full codons
    aa = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    # stop at first stop codon if present (avoid translating past it)
    if '*' in aa:
        aa = aa[:aa.index('*')]
    return "".join(aa)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out-protein-fasta", required=True)
    parser.add_argument("--out-nt-alignment", required=True)
    parser.add_argument("--mafft-args", default="--auto",
                         help="Extra args passed to mafft for the protein alignment")
    args = parser.parse_args()

    nt_seqs = parse_fasta(args.fasta)
    protein_seqs = {}
    skipped = []
    for header, seq in nt_seqs.items():
        prot = translate(seq)
        if len(prot) < 10:
            skipped.append(header)
            continue
        protein_seqs[header] = prot

    print(f"Translated {len(protein_seqs)}/{len(nt_seqs)} sequences "
          f"({len(skipped)} skipped, too short after stop-codon trim)")

    with open(args.out_protein_fasta, "w") as f:
        for h, p in protein_seqs.items():
            f.write(f">{h}\n{p}\n")

    print(f"Running MAFFT on protein alignment...")
    with open(args.out_protein_fasta.replace(".fasta", "_aligned.fasta"), "w") as out:
        result = subprocess.run(
            ["mafft"] + args.mafft_args.split() + [args.out_protein_fasta],
            stdout=out, stderr=subprocess.PIPE, text=True
        )
    if result.returncode != 0:
        print("MAFFT failed:", result.stderr, file=sys.stderr)
        sys.exit(1)

    protein_aligned = parse_fasta(args.out_protein_fasta.replace(".fasta", "_aligned.fasta"))
    aln_len = len(next(iter(protein_aligned.values())))
    print(f"Protein alignment: {len(protein_aligned)} sequences x {aln_len} columns")

    # Map back to codon-aligned nucleotide alignment.
    with open(args.out_nt_alignment, "w") as out:
        for header, prot_aligned in protein_aligned.items():
            nt_seq = nt_seqs[header]
            codon_aligned = []
            nt_pos = 0
            for aa in prot_aligned:
                if aa == "-":
                    codon_aligned.append("---")
                else:
                    codon_aligned.append(nt_seq[nt_pos:nt_pos + 3])
                    nt_pos += 3
            out.write(f">{header}\n{''.join(codon_aligned)}\n")

    print(f"Codon-aligned nucleotide alignment written to {args.out_nt_alignment}")


if __name__ == "__main__":
    main()
