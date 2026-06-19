"""
Prepare clean training data for BGC-SetNet.
Applies three critical fixes:
  1. Require >= 3 genes hitting same MIBiG reference for compound family assignment
  2. Deduplicate near-identical BGC variants within each family (Jaccard > 0.8)
  3. Assign product_class from probability columns where missing
"""
import pandas as pd
import numpy as np
from collections import defaultdict


def assign_product_class(atlas):
    """Fix 3: Fill missing product_class using argmax of probability columns."""
    prob_cols = ['Alkaloid', 'NRP', 'Other', 'Polyketide', 'RiPP', 'Saccharide', 'Terpene']
    missing = atlas['product_class'].isna()
    atlas.loc[missing, 'product_class'] = atlas.loc[missing, prob_cols].idxmax(axis=1)
    return atlas


def assign_compound_families(genes, atlas, min_genes=3):
    """Fix 1: Only assign compound family if >= min_genes from same BGC hit same MIBiG ref."""
    hits = genes[(genes['mibig_id'] != 'No Match') & (genes['bgc_id'].isin(set(atlas['bgc_id'])))].copy()
    hits['mibig_ref'] = hits['mibig_id'].str.extract(r'(BGC\d+)')[0]

    # Count genes per (BGC, MIBiG reference) pair
    pair_counts = hits.groupby(['bgc_id', 'mibig_ref']).size().reset_index(name='n_gene_hits')

    # Keep only pairs with >= min_genes hits
    strong_pairs = pair_counts[pair_counts['n_gene_hits'] >= min_genes]

    # For each BGC, pick the MIBiG ref with the most gene hits (best match)
    best_ref = strong_pairs.sort_values('n_gene_hits', ascending=False).drop_duplicates('bgc_id')
    best_ref = best_ref[['bgc_id', 'mibig_ref', 'n_gene_hits']]
    best_ref.columns = ['bgc_id', 'compound_family', 'n_supporting_genes']

    return best_ref


def pfam_jaccard(pfam_set_a, pfam_set_b):
    """Jaccard similarity between two Pfam domain sets."""
    if not pfam_set_a and not pfam_set_b:
        return 1.0
    intersection = len(pfam_set_a & pfam_set_b)
    union = len(pfam_set_a | pfam_set_b)
    return intersection / union if union > 0 else 0.0


def deduplicate_family(family_bgcs, atlas, threshold=0.8):
    """
    Fix 2: Within a compound family, collapse near-identical BGC variants.
    Returns representative BGC IDs and marks trivial vs informative pairs.

    Two BGCs with Pfam Jaccard > threshold are considered 'same variant' (trivial copy).
    Keep one representative per variant group.
    """
    if len(family_bgcs) <= 1:
        return family_bgcs, []

    bgc_pfams = {}
    for bgc_id in family_bgcs:
        row = atlas[atlas['bgc_id'] == bgc_id]
        if row.empty:
            continue
        pfam_str = row.iloc[0]['pfam_ids']
        bgc_pfams[bgc_id] = set(str(pfam_str).split(';')) if pd.notna(pfam_str) else set()

    # Greedy clustering: assign each BGC to first cluster it's similar enough to
    clusters = []  # list of lists
    for bgc_id, pfams in bgc_pfams.items():
        placed = False
        for cluster in clusters:
            rep_pfams = bgc_pfams[cluster[0]]
            if pfam_jaccard(pfams, rep_pfams) > threshold:
                cluster.append(bgc_id)
                placed = True
                break
        if not placed:
            clusters.append([bgc_id])

    # Representatives: one per cluster (pick the one with most proteins)
    representatives = []
    for cluster in clusters:
        cluster_rows = atlas[atlas['bgc_id'].isin(cluster)]
        if cluster_rows.empty:
            representatives.append(cluster[0])
        else:
            best = cluster_rows.loc[cluster_rows['num_proteins'].idxmax(), 'bgc_id']
            representatives.append(best)

    # Informative pairs: representatives from DIFFERENT clusters (architecturally distinct)
    informative_pairs = []
    for i in range(len(representatives)):
        for j in range(i + 1, len(representatives)):
            informative_pairs.append((representatives[i], representatives[j]))

    return representatives, informative_pairs


def add_gene_position(genes, atlas):
    """Fix 4: Add relative position (0-1) of each gene within its BGC."""
    atlas_coords = atlas.set_index('bgc_id')[['nucl_start', 'nucl_end']].to_dict('index')

    positions = []
    for _, row in genes.iterrows():
        bgc_id = row['bgc_id']
        if bgc_id in atlas_coords:
            bgc_start = atlas_coords[bgc_id]['nucl_start']
            bgc_end = atlas_coords[bgc_id]['nucl_end']
            bgc_length = bgc_end - bgc_start
            if bgc_length > 0:
                # Extract gene position from gene_id (format: scaffold_scaffold_GENENUM)
                # Approximate: use gene index within BGC / total genes
                positions.append(np.nan)  # Will compute from ordering below
            else:
                positions.append(0.5)
        else:
            positions.append(0.5)

    # Better approach: rank genes within each BGC by their order in the file
    # (genes are listed in genomic order in the CSV)
    genes = genes.copy()
    genes['gene_rank'] = genes.groupby('bgc_id').cumcount()
    genes['gene_total'] = genes.groupby('bgc_id')['bgc_id'].transform('count')
    genes['relative_position'] = genes['gene_rank'] / genes['gene_total'].clip(lower=1)

    return genes


def main():
    print("Loading data...")
    atlas = pd.read_csv('community_atlas.csv')
    discovery = pd.read_csv('discovery_zone.csv')
    genes = pd.read_csv('master_gene_annotations.csv')

    print(f"Starting: {len(atlas)} atlas BGCs, {len(discovery)} discovery BGCs, {len(genes)} genes")
    print()

    # === FIX 3: Fill missing product_class ===
    missing_before = atlas['product_class'].isna().sum()
    atlas = assign_product_class(atlas)
    print(f"Fix 3 - Product class: filled {missing_before} missing values")
    print(f"  Distribution: {atlas['product_class'].value_counts().to_dict()}")
    print()

    # === FIX 1: Strict compound family assignment ===
    family_assignments = assign_compound_families(genes, atlas, min_genes=3)
    print(f"Fix 1 - Compound families (>= 3 genes required):")
    print(f"  BGCs with reliable family label: {len(family_assignments)} / {len(atlas)}")
    print(f"  Unique compound families: {family_assignments['compound_family'].nunique()}")
    print(f"  Mean supporting genes: {family_assignments['n_supporting_genes'].mean():.1f}")
    print()

    # Merge family labels into atlas
    atlas = atlas.merge(family_assignments, on='bgc_id', how='left')
    labeled = atlas[atlas['compound_family'].notna()]
    unlabeled = atlas[atlas['compound_family'].isna()]
    print(f"  Labeled (trainable with contrastive loss): {len(labeled)}")
    print(f"  Unlabeled (use for pre-training only): {len(unlabeled)}")
    print()

    # === FIX 2: Deduplicate within families ===
    print("Fix 2 - Deduplicating near-identical variants (Jaccard > 0.8)...")
    families = labeled.groupby('compound_family')['bgc_id'].apply(list).to_dict()

    all_representatives = []
    all_informative_pairs = []
    total_before = 0
    total_after = 0

    for fam_id, bgc_list in families.items():
        total_before += len(bgc_list)
        reps, pairs = deduplicate_family(bgc_list, atlas, threshold=0.8)
        total_after += len(reps)
        for r in reps:
            all_representatives.append({'bgc_id': r, 'compound_family': fam_id})
        for p in pairs:
            all_informative_pairs.append({
                'anchor': p[0], 'positive': p[1], 'compound_family': fam_id
            })

    print(f"  Before dedup: {total_before} BGCs across {len(families)} families")
    print(f"  After dedup: {total_after} representative BGCs")
    print(f"  Informative positive pairs: {len(all_informative_pairs)}")
    print(f"  Reduction: {(1 - total_after/total_before)*100:.1f}%")
    print()

    # === FIX 4: Add gene position ===
    print("Fix 4 - Adding relative gene positions...")
    genes = add_gene_position(genes, atlas)
    print(f"  Added relative_position column (0.0 = start of cluster, 1.0 = end)")
    print()

    # === SAVE OUTPUTS ===
    # Training-ready datasets
    reps_df = pd.DataFrame(all_representatives)
    pairs_df = pd.DataFrame(all_informative_pairs)

    atlas.to_csv('community_atlas.csv', index=False)
    genes.to_csv('master_gene_annotations.csv', index=False)
    reps_df.to_csv('training_representatives.csv', index=False)
    pairs_df.to_csv('training_positive_pairs.csv', index=False)

    print("=== SAVED FILES ===")
    print(f"  community_atlas.csv - updated with compound_family + filled product_class")
    print(f"  master_gene_annotations.csv - updated with relative_position")
    print(f"  training_representatives.csv - {len(reps_df)} deduplicated BGCs for training")
    print(f"  training_positive_pairs.csv - {len(pairs_df)} informative pairs (architecturally distinct)")
    print()

    # === SUMMARY STATS ===
    trainable_families = reps_df.groupby('compound_family').size()
    print("=== FINAL TRAINING DATA SUMMARY ===")
    print(f"  Total representative BGCs: {len(reps_df)}")
    print(f"  Compound families: {len(trainable_families)}")
    print(f"  Families with >= 2 reps (can form pairs): {(trainable_families >= 2).sum()}")
    print(f"  Families with >= 5 reps: {(trainable_families >= 5).sum()}")
    print(f"  Informative positive pairs: {len(pairs_df)}")
    print(f"  Unlabeled BGCs (for pre-training): {len(unlabeled) + len(discovery)}")
    print(f"  Discovery zone (for inference): {len(discovery)}")


if __name__ == '__main__':
    main()
