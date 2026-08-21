# Bioinformatics Data Pipeline

This folder contains the lightweight, source-controlled pieces of the BGC data-preparation workflow used upstream of the retrieval benchmark. Large raw and generated artifacts are intentionally excluded from Git.

## What Is Tracked

- `scripts/run_deepbgc_worker.sh`: one-genome DeepBGC worker suitable for GNU parallel or a scheduler wrapper.
- `scripts/build_atlas_tables.py`: merges gene-level BLASTp hits with DeepBGC BGC metadata and emits `master_gene_annotations.csv`, `community_atlas.csv`, and `discovery_zone.csv`.
- `split.py`: the original notebook-era atlas split script retained for provenance. Prefer `scripts/build_atlas_tables.py` for new runs because it exposes thresholds as CLI arguments.
- `antismash/assembly_ids.txt` and `antismash/test_ids.txt`: small assembly ID manifests used during antiSMASH/DeepBGC checks.

## Generated Inputs Not Tracked

The following files are required to rerun the full data build but are too large or too derived for this Git repository:

- genome FASTA/GenBank files
- antiSMASH and DeepBGC result directories
- BLAST databases and raw BLAST result tables
- `all_genes.csv` and `all_genes.fasta`
- `all_bgcs_combined.csv`
- `master_gene_annotations.csv`
- `community_atlas.csv`
- `discovery_zone.csv`
- `esm2_embeddings.h5`

The final five-seed model checkpoints are released separately on Hugging Face: https://huggingface.co/whiteh4t/bgc-setnet

## Pipeline Sketch

1. Deduplicate genome assembly IDs with `deduplicate_genomes.py`.
2. Run DeepBGC per genome, for example with `bioinformatics/scripts/run_deepbgc_worker.sh`.
3. Export all detected BGC gene records into `all_genes.csv` and protein sequences into `all_genes.fasta`.
4. BLAST proteins against MIBiG proteins to create `blast_raw_results.tsv` with columns `gene_id`, `mibig_id`, `pident`, `qcovs`, `evalue`.
5. Build atlas tables:

```bash
python bioinformatics/scripts/build_atlas_tables.py \
  --genes all_genes.csv \
  --blast blast_raw_results.tsv \
  --bgcs all_bgcs_combined.csv \
  --output-dir .
```

6. Prepare stricter downstream training/evaluation labels with `prepare_training_data.py` where needed.
7. Filter proteins and compute frozen ESM-2 embeddings with `filter_proteins.py` and `esm2_inference.py`.
8. Use the final-paper source/model artifacts released on Hugging Face for the five-seed retrieval experiments.

## Default Thresholds

`build_atlas_tables.py` uses the historical preprocessing thresholds:

- minimum BLAST query coverage: `50.0`
- known/community identity threshold: `40.0`
- minimum DeepBGC score: `0.7`

These thresholds reproduce the atlas/discovery-zone construction logic, not the final grouped retrieval split by themselves. The final split and model artifacts are documented in the GitHub root README, `CHECKPOINTS.md`, and the Hugging Face model card.
