"""
Filter all_genes.fasta to only proteins present in master_gene_annotations.csv.
Produces filtered_proteins.fasta ready for ESM-2 inference.

Input:  ../bioinformatics/all_genes.fasta (~532K proteins)
Output: data/filtered_proteins.fasta (316,257 proteins)
"""
import os
import pandas as pd

FASTA_IN = "../bioinformatics/all_genes.fasta"
CSV_PATH = "master_gene_annotations.csv"
OUT_DIR = "data"
FASTA_OUT = os.path.join(OUT_DIR, "filtered_proteins.fasta")


def load_gene_ids(csv_path):
    df = pd.read_csv(csv_path, usecols=["gene_id"])
    return set(df["gene_id"].values)


def filter_fasta(fasta_in, fasta_out, keep_ids):
    kept = 0
    skipped = 0
    seen = set()
    with open(fasta_in, "r") as fin, open(fasta_out, "w") as fout:
        writing = False
        for line in fin:
            if line.startswith(">"):
                gene_id = line[1:].strip()
                if gene_id in keep_ids and gene_id not in seen:
                    writing = True
                    kept += 1
                    seen.add(gene_id)
                    fout.write(line)
                else:
                    writing = False
                    skipped += 1
            elif writing:
                fout.write(line)
    return kept, skipped


def main():
    print("Loading gene IDs from master_gene_annotations.csv...")
    keep_ids = load_gene_ids(CSV_PATH)
    print(f"  Gene IDs to keep: {len(keep_ids):,}")

    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"\nFiltering {FASTA_IN} → {FASTA_OUT}...")
    kept, skipped = filter_fasta(FASTA_IN, FASTA_OUT, keep_ids)

    print(f"\nDone.")
    print(f"  Kept:    {kept:,} proteins")
    print(f"  Skipped: {skipped:,} proteins")
    print(f"  Output:  {FASTA_OUT}")

    if kept != len(keep_ids):
        missing = len(keep_ids) - kept
        print(f"\n  WARNING: {missing:,} gene IDs from CSV not found in FASTA!")


if __name__ == "__main__":
    main()
