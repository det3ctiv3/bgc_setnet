import pandas as pd

print("Loading data...")
genes_df = pd.read_csv('all_genes.csv')
blast_cols = ['gene_id', 'mibig_id', 'pident', 'qcovs', 'evalue']
blast_df = pd.read_csv('blast_raw_results.tsv', sep='\t', names=blast_cols)
bgc_metadata = pd.read_csv('all_bgcs_combined.csv')

valid_blast = blast_df[blast_df['qcovs'] >= 50.0]

print('Merging datasets...')
merged_genes = pd.merge(genes_df, valid_blast, on='gene_id', how='left')
merged_genes['pident'] = merged_genes['pident'].fillna(0)
merged_genes['mibig_id'] = merged_genes['mibig_id'].fillna("No Match")
merged_genes.to_csv('master_gene_annotations.csv', index=False)

print('Calculating BGC novelty scores...')
hits_only = merged_genes[merged_genes['pident'] > 0]

bgc_novelty = hits_only.groupby(['genome_id', 'bgc_id'])['pident'].mean().reset_index()
bgc_novelty.rename(columns={'pident': 'avg_mibig_identity'}, inplace=True)

all_bgcs = genes_df[['genome_id', 'bgc_id']].drop_duplicates()
bgc_novelty = pd.merge(all_bgcs, bgc_novelty, on=['genome_id', 'bgc_id'], how='left')
bgc_novelty['avg_mibig_identity'] = bgc_novelty['avg_mibig_identity'].fillna(0)

print('Integrating DeepBGC Confidence Scores...')
if 'bgc_id' not in bgc_metadata.columns and 'bgc_candidate_id' in bgc_metadata.columns:
    bgc_metadata.rename(columns={'bgc_candidate_id': 'bgc_id'}, inplace=True)

final_bgc_table = pd.merge(bgc_novelty, bgc_metadata, on=['genome_id', 'bgc_id'], how='inner')
score_col = 'deepbgc_score' if 'deepbgc_score' in final_bgc_table.columns else 'score'

community_atlas = final_bgc_table[
    (final_bgc_table['avg_mibig_identity'] >= 40.0) &
    (final_bgc_table[score_col] >= 0.7)
]

discovery_zone = final_bgc_table[
    (final_bgc_table['avg_mibig_identity'] < 40.0) &
    (final_bgc_table[score_col] >= 0.7)
]

community_atlas.to_csv('community_atlas.csv', index=False)
discovery_zone.to_csv('discovery_zone.csv', index=False)

print(f"Total BGCs processed: {len(final_bgc_table)}")
print(f"Community Atlas (Known Pathways): {len(community_atlas)} BGCs")
print(f"Discovery Zone (True Novel Candidates): {len(discovery_zone)} BGCs")
