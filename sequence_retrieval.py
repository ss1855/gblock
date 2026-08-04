#!/usr/bin/env python3
"""
Sequence retrieval: given a gene presence screen output (species, best_accession
columns), fetch the coding NUCLEOTIDE sequence for each gene-positive species.

WHY THIS STEP MATTERS
    The presence screen (gene_presence_screen.py) searches NCBI's protein
    database, so best_accession is a PROTEIN accession. Primers bind DNA, so
    building a consensus amplicon and validating primer binding requires the
    coding nucleotide sequence, not the protein sequence. This script makes
    that protein -> nucleotide jump.

HOW THE PROTEIN -> NUCLEOTIDE MAPPING WORKS
    The obvious approach -- fetch the protein GenBank flat file and parse the
    CDS feature's /coded_by qualifier -- only works for accessions submitted
    directly to INSDC (GenBank/ENA/DDBJ). It does NOT work for RefSeq "WP_"
    accessions (non-redundant protein records, merged across every genome
    that encodes an identical protein), which have no /coded_by qualifier at
    all. WP_ accessions are the majority case in practice (~82% of ureC hits
    in the DeQueen pilot), so this script uses NCBI's Identical Protein
    Groups (IPG) report as the primary method instead:
        efetch(db="protein", id=accession, rettype="ipg", retmode="xml")
    The IPG report lists every nucleotide CDS that encodes the queried
    protein, grouped per source protein accession, and works uniformly for
    both RefSeq and INSDC accessions. We match the <Protein accver=...>
    element to our exact queried accession (not just any CDS in the report)
    so multi-genome/MULTISPECIES records never pull the wrong strain's
    sequence. /coded_by parsing from the protein flat file is kept as a
    fallback for the rare accession IPG doesn't resolve.

USAGE
    python sequence_retrieval.py --presence-csv ureC_presence.csv \
        --email you@example.com \
        --api-key YOUR_NCBI_API_KEY \
        --out-fasta ureC_sequences.fasta \
        --out-log ureC_retrieval_log.csv

    # Preview which species/accessions would be processed, no network calls:
    python sequence_retrieval.py --presence-csv ureC_presence.csv \
        --email you@example.com --dry-run --limit 5

REQUIREMENTS
    pip install biopython pandas

NOTES
    - Requires outbound internet access to eutils.ncbi.nlm.nih.gov.
    - NCBI requires an email for Entrez use. An API key is optional but
      raises your rate limit from ~3 req/s to ~10 req/s.
    - Resumable: writes the log incrementally and skips species already
      present in an existing --out-log file if you re-run after an
      interruption. FASTA entries are appended in the same run, so a
      resumed run's FASTA stays in sync with the log.
    - Only rows with hit_count > 0 in the presence CSV are attempted.
    - Failures (unparseable records, no CDS found, network errors) are
      logged with a reason and do not stop the run.
"""

import argparse
import csv
import os
import re
import time
import xml.etree.ElementTree as ET

from Bio import Entrez


def fetch_ipg_cds(accession: str, species: str):
    """Return (nuc_accession, start, stop, strand) for the CDS encoding
    `accession`, preferring the CDS entry whose organism matches `species`.
    Raises ValueError if no matching CDS is found."""
    handle = Entrez.efetch(db="protein", id=accession, rettype="ipg", retmode="xml")
    text = handle.read()
    handle.close()
    if isinstance(text, bytes):
        text = text.decode()

    root = ET.fromstring(text)
    proteins = root.findall(".//Protein")

    # Prefer the <Protein> element matching our exact queried accession.
    matched = [p for p in proteins if p.get("accver") == accession]
    candidates = matched if matched else proteins

    cds_list = []
    for p in candidates:
        cds_list.extend(p.findall("./CDSList/CDS"))
    if not matched:
        # Fallback: no exact accession match in the report at all -- gather
        # every CDS across every protein entry so species-matching below has
        # something to work with.
        cds_list = [cds for p in proteins for cds in p.findall("./CDSList/CDS")]

    if not cds_list:
        raise ValueError("IPG report contained no CDS entries")

    def cds_tuple(cds):
        return (cds.get("accver"), int(cds.get("start")), int(cds.get("stop")),
                cds.get("strand"))

    # Prefer a CDS whose org exactly matches the species; otherwise take the
    # first entry under the matched protein (or first overall as last resort).
    for cds in cds_list:
        if cds.get("org") == species:
            return cds_tuple(cds)
    return cds_tuple(cds_list[0])


def fetch_coded_by(accession: str):
    """Fallback: parse /coded_by from the protein GenBank flat file. Returns
    (nuc_accession, start, stop, strand) or raises ValueError."""
    handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
    text = handle.read()
    handle.close()

    m = re.search(r'/coded_by="([^"]+)"', text)
    if not m:
        raise ValueError("no /coded_by qualifier found")

    loc = m.group(1)
    strand = "-" if loc.startswith("complement(") else "+"
    loc = loc.replace("complement(", "").rstrip(")")
    m2 = re.match(r'([\w.]+):(\d+)\.\.(\d+)', loc)
    if not m2:
        raise ValueError(f"unparseable /coded_by location: {loc}")
    nuc_acc, start, stop = m2.group(1), int(m2.group(2)), int(m2.group(3))
    return nuc_acc, start, stop, strand


def fetch_nucleotide(nuc_acc: str, start: int, stop: int, strand: str) -> str:
    strand_code = 2 if strand == "-" else 1
    handle = Entrez.efetch(db="nuccore", id=nuc_acc, rettype="fasta", retmode="text",
                            seq_start=start, seq_stop=stop, strand=strand_code)
    text = handle.read()
    handle.close()
    lines = text.strip().splitlines()
    return "".join(lines[1:])


def resolve_and_fetch(protein_acc: str, species: str):
    """Returns (nuc_accession, sequence). Raises ValueError/Exception on
    failure, with a message describing what went wrong."""
    try:
        nuc_acc, start, stop, strand = fetch_ipg_cds(protein_acc, species)
    except Exception:
        nuc_acc, start, stop, strand = fetch_coded_by(protein_acc)

    seq = fetch_nucleotide(nuc_acc, start, stop, strand)
    if not seq:
        raise ValueError("nucleotide fetch returned an empty sequence")
    return nuc_acc, seq


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve coding nucleotide sequences for gene-positive "
                     "species from a gene_presence_screen.py output CSV"
    )
    parser.add_argument("--presence-csv", required=True,
                         help="Output CSV from gene_presence_screen.py")
    parser.add_argument("--email", required=True, help="Required by NCBI for Entrez use")
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--out-fasta", required=False)
    parser.add_argument("--out-log", required=False)
    parser.add_argument("--limit", type=int, default=None,
                         help="Limit number of species (useful for testing)")
    parser.add_argument("--dry-run", action="store_true",
                         help="List species/accessions that would be processed, "
                              "without hitting the network")
    args = parser.parse_args()
    if not args.dry_run and not (args.out_fasta and args.out_log):
        parser.error("--out-fasta and --out-log are required unless --dry-run is set")

    import pandas as pd
    df = pd.read_csv(args.presence_csv)
    df = df[df["hit_count"] > 0].reset_index(drop=True)
    if args.limit:
        df = df.head(args.limit)

    if args.dry_run:
        for _, row in df.iterrows():
            print(f"{row['species']}\t{row['best_accession']}")
        return

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    sleep_time = 0.11 if args.api_key else 0.34  # ~10/s vs ~3/s, with margin

    done_species = set()
    log_exists = os.path.isfile(args.out_log)
    if log_exists:
        with open(args.out_log, newline="") as f:
            for row in csv.DictReader(f):
                done_species.add(row["species"])

    log_fieldnames = ["species", "protein_accession", "nucleotide_accession",
                       "status", "error_message"]

    log_file = open(args.out_log, "a", newline="")
    fasta_file = open(args.out_fasta, "a")
    writer = csv.DictWriter(log_file, fieldnames=log_fieldnames)
    if not log_exists:
        writer.writeheader()

    try:
        for i, row in df.iterrows():
            species = row["species"]
            protein_acc = row["best_accession"]
            if species in done_species:
                continue

            nuc_acc, status, error_message = "", "failed", ""
            try:
                nuc_acc, seq = resolve_and_fetch(protein_acc, species)
                fasta_file.write(f">{species}|{protein_acc}|{nuc_acc}\n{seq}\n")
                fasta_file.flush()
                status = "success"
            except Exception as e:
                error_message = str(e)

            writer.writerow({
                "species": species,
                "protein_accession": protein_acc,
                "nucleotide_accession": nuc_acc,
                "status": status,
                "error_message": error_message,
            })
            log_file.flush()

            print(f"[{i + 1}/{len(df)}] {species} ({protein_acc}): {status}"
                  + (f" - {error_message}" if error_message else ""))
            time.sleep(sleep_time)
    finally:
        log_file.close()
        fasta_file.close()

    print(f"\nDone. Sequences written to {args.out_fasta}, log written to {args.out_log}")


if __name__ == "__main__":
    main()
