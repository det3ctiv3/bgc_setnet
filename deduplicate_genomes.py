"""
Deduplicate genomes: when both GCA_ and GCF_ exist for the same assembly,
keep only GCF_ (RefSeq = higher curation). Keep GCA_ only when no GCF_ exists.
"""
import pandas as pd

def get_base_id(genome_id):
    return genome_id.replace("GCA_", "").replace("GCF_", "")

def build_keep_set(genome_ids):
    base_to_genomes = {}
    for gid in genome_ids:
        base = get_base_id(gid)
        base_to_genomes.setdefault(base, []).append(gid)

    keep = set()
    for base, gids in base_to_genomes.items():
        gcf = [g for g in gids if g.startswith("GCF_")]
        gca = [g for g in gids if g.startswith("GCA_")]
        if gcf:
            keep.add(gcf[0])
        elif gca:
            keep.add(gca[0])
    return keep

# Load all three CSVs
print("Loading CSVs...")
atlas = pd.read_csv("community_atlas.csv")
discovery = pd.read_csv("discovery_zone.csv")
genes = pd.read_csv("master_gene_annotations.csv")

# Get all genome IDs across datasets
all_genomes = set(atlas["genome_id"]).union(discovery["genome_id"]).union(genes["genome_id"])
print(f"Total genome entries before dedup: {len(all_genomes)}")

# Build keep set
keep = build_keep_set(all_genomes)
print(f"Unique genomes after dedup: {len(keep)}")
print(f"Removed duplicates: {len(all_genomes) - len(keep)}")

# Filter
atlas_clean = atlas[atlas["genome_id"].isin(keep)]
discovery_clean = discovery[discovery["genome_id"].isin(keep)]
genes_clean = genes[genes["genome_id"].isin(keep)]

print(f"\ncommunity_atlas: {len(atlas)} -> {len(atlas_clean)} BGCs")
print(f"discovery_zone: {len(discovery)} -> {len(discovery_clean)} BGCs")
print(f"master_gene_annotations: {len(genes)} -> {len(genes_clean)} genes")

# Save cleaned versions
atlas_clean.to_csv("community_atlas.csv", index=False)
discovery_clean.to_csv("discovery_zone.csv", index=False)
genes_clean.to_csv("master_gene_annotations.csv", index=False)

print("\nDone. Original files overwritten with deduplicated versions.")
