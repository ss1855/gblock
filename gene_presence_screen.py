#!/usr/bin/env python3
"""
Gene presence screen: determine which DeQueen bacterial species have a given
target gene, by querying NCBI's protein database via Entrez.

WHY THIS STEP MATTERS
    ureC, amoA, nirK, nirS, and nosZ are functional marker genes, not
    universal genes like 16S. Most of the 1584 DeQueen species will NOT
    carry any given target. Building a consensus amplicon from the full
    species list (as if every organism had the gene) would misrepresent
    presence/absence. This script narrows each gene down to its actual
    gene-positive subset first.

USAGE
    python gene_presence_screen.py --gene ureC \
        --tax-table tax_table_df.csv \
        --email you@example.com \
        --api-key YOUR_NCBI_API_KEY \
        --out ureC_presence.csv

    # Preview the search queries without hitting the network:
    python gene_presence_screen.py --gene ureC --tax-table tax_table_df.csv \
        --email you@example.com --dry-run --limit 5

REQUIREMENTS
    pip install biopython pandas

NOTES
    - Requires outbound internet access to eutils.ncbi.nlm.nih.gov. This was
      written and tested in a sandbox that cannot reach NCBI (confirmed
      403 from the egress proxy) -- run this on your AWS/SageMaker
      environment or any machine with normal internet access.
    - NCBI requires an email for Entrez use. An API key is optional but
      raises your rate limit from ~3 req/s to ~10 req/s:
      https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
    - Resumable: writes results incrementally and skips species already
      present in an existing --out file if you re-run after an interruption.
    - This screens by NCBI protein-record title text, which is a recall
      step, not a definitive presence/absence call. Hits should be spot
      checked before being treated as ground truth, and a hit_count of 0
      does not prove the gene is absent, only that no annotated protein
      record with that name was found for that organism.
"""

import argparse
import csv
import os
import time

from Bio import Entrez

# Gene name aliases, widen search recall since NCBI records use inconsistent
# naming across submitters (gene symbol vs. full protein name).
GENE_ALIASES = {
    "ureC": ["ureC", "urease subunit alpha", "urease alpha subunit"],
    "amoA_AOB": ["amoA", "ammonia monooxygenase subunit A"],
    "amoA_AOA": ["amoA", "ammonia monooxygenase subunit A"],
    "nirK": ["nirK", "copper-containing nitrite reductase",
             "nitrite reductase, copper-containing"],
    "nirS": ["nirS", "cytochrome cd1 nitrite reductase",
             "nitrite reductase, cytochrome cd1"],
    "nosZ_cladeI": ["nosZ", "nitrous oxide reductase"],
    "nosZ_cladeII": ["nosZ", "nitrous oxide reductase"],
}


def build_query(gene: str, organism: str) -> str:
    aliases = GENE_ALIASES.get(gene, [gene])
    alias_terms = " OR ".join(f'"{a}"[Title]' for a in aliases)
    return f'({alias_terms}) AND "{organism}"[Organism]'


def search_gene(gene: str, organism: str, retmax: int = 5):
    query = build_query(gene, organism)
    handle = Entrez.esearch(db="protein", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record


def fetch_summary(ids):
    if not ids:
        return []
    handle = Entrez.esummary(db="protein", id=",".join(ids))
    records = Entrez.read(handle)
    handle.close()
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Screen DeQueen species for target gene presence via NCBI"
    )
    parser.add_argument("--gene", required=True, choices=list(GENE_ALIASES.keys()))
    parser.add_argument("--tax-table", required=True)
    parser.add_argument("--email", required=True, help="Required by NCBI for Entrez use")
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--out", required=False)
    parser.add_argument("--limit", type=int, default=None,
                         help="Limit number of species (useful for testing)")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print queries without hitting the network")
    args = parser.parse_args()
    if not args.dry_run and not args.out:
        parser.error("--out is required unless --dry-run is set")

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key
    sleep_time = 0.11 if args.api_key else 0.34  # ~10/s vs ~3/s, with margin

    import pandas as pd
    df = pd.read_csv(args.tax_table)
    df = df[df["Domain"] == "Bacteria"].reset_index(drop=True)
    if args.limit:
        df = df.head(args.limit)

    if args.dry_run:
        for _, row in df.iterrows():
            print(build_query(args.gene, row["Species"]))
        return

    done_species = set()
    file_exists = os.path.isfile(args.out)
    if file_exists:
        with open(args.out, newline="") as f:
            for row in csv.DictReader(f):
                done_species.add(row["species"])

    fieldnames = ["species", "genus", "phylum", "gene_target", "hit_count",
                  "best_accession", "best_title", "query"]

    with open(args.out, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()

        for i, row in df.iterrows():
            species = row["Species"]
            if species in done_species or species == "Unknown":
                continue

            query = build_query(args.gene, species)
            count, best_accession, best_title = -1, "", ""
            try:
                result = search_gene(args.gene, species)
                count = int(result["Count"])
                ids = result["IdList"]
                if ids:
                    summaries = fetch_summary(ids[:1])
                    if summaries:
                        best_accession = summaries[0].get("AccessionVersion", "")
                        best_title = summaries[0].get("Title", "")
            except Exception as e:
                best_accession = f"ERROR: {e}"

            writer.writerow({
                "species": species,
                "genus": row["Genus"],
                "phylum": row["Phylum"],
                "gene_target": args.gene,
                "hit_count": count,
                "best_accession": best_accession,
                "best_title": best_title,
                "query": query,
            })
            f.flush()

            print(f"[{i + 1}/{len(df)}] {species}: {count} hits")
            time.sleep(sleep_time)

    print(f"\nDone. Results written to {args.out}")


if __name__ == "__main__":
    main()
