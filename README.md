# BGC-SetNet Retrieval Benchmark

Analysis code, result artifacts, and manuscript source for:

> Sequence-Derived Representations versus Pfam-Domain Content for Biosynthetic Gene Cluster Retrieval

This repository reports a group-disjoint, five-seed evaluation of biosynthetic gene cluster (BGC) retrieval in a *Streptomyces griseus* species-family atlas. The central result is negative but informative: explicit Pfam-domain content remains the strongest signal on this benchmark.

## Final Results

The effective inference unit is the held-out family (n = 16). Values below are five-seed family means.

| Method | Recall@50 | MRR | MAP | nDCG@50 |
| --- | ---: | ---: | ---: | ---: |
| Raw ESM mean | 0.7946 | 0.2550 | 0.7251 | 0.8078 |
| BGC-SetNet + Pfam | 0.8472 | 0.2786 | 0.7771 | 0.8502 |
| Pfam Jaccard | 0.8788 | 0.3071 | 0.8480 | 0.9042 |
| ESM + BGC-SetNet + Pfam | 0.8769 | 0.3096 | 0.8503 | 0.9058 |
| Weighted Pfam Jaccard | 0.8789 | 0.3069 | 0.8477 | 0.9040 |

Weighted Pfam exceeds unweighted Pfam by only 0.00003 Recall@50 (p = 0.50). The ESM plus Pfam-augmented BGC-SetNet ensemble is lower by 0.00198 (p = 0.625). Neither difference supports an improvement claim.

## Release Scope

This release contains:

- `bioinformatics/`: curated data-preparation pipeline notes, DeepBGC worker script, atlas-table builder, and small assembly ID manifests.
- `deduplicate_genomes.py`, `filter_proteins.py`, `esm2_inference.py`, `prepare_training_data.py`: preprocessing utilities used around the atlas/model workflow.
- `build_final_results.py`: aggregation, bootstrap intervals, and exact paired tests.
- `plot_final_results.py`: deterministic regeneration of the four manuscript figures.
- `results/dgx_final/`: compact canonical aggregate and per-seed family results.
- `paper/main.tex`: canonical manuscript source.
- `paper/fig_final_*.pdf`: figures generated from the canonical result tables.
- `CHECKPOINTS.md`: expected source checkpoint hashes used to verify the public Hugging Face release.

Large raw/generated assets are intentionally excluded from Git: genome FASTA/GenBank files, antiSMASH/DeepBGC output directories, BLAST databases/results, full atlas CSVs, ESM-2 embedding HDF5 files, and local legacy checkpoints.

## Reproduce Analysis Figures

Use Python 3.10 or newer.

~~~bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python build_final_results.py
python plot_final_results.py
~~~

The committed results/dgx_final/ tables are sufficient for plot_final_results.py. Rebuilding the aggregate with build_final_results.py additionally requires the retained v5/v6 per-seed evaluation exports.

## Reproduce Data Preparation

The curated upstream data pipeline is documented in `bioinformatics/README.md`. In short, genome assemblies are processed with DeepBGC, protein hits are compared against MIBiG with BLASTp, atlas/discovery tables are generated with `bioinformatics/scripts/build_atlas_tables.py`, and frozen ESM-2 gene embeddings are computed with `esm2_inference.py`. The large generated inputs and outputs are not stored in Git; model checkpoints and per-seed evaluation artifacts are on Hugging Face: https://huggingface.co/whiteh4t/bgc-setnet.

Build the manuscript from the canonical paper directory:

~~~bash
cd paper
tectonic main.tex
~~~

## Reproducibility Boundary

The exact v5 Pfam-augmented BGC-SetNet checkpoints, v6 weighted-Pfam checkpoints, training-only Pfam vocabularies, per-seed evaluation artifacts, source snapshot, and model card are publicly released on Hugging Face: https://huggingface.co/whiteh4t/bgc-setnet. The available July checkpoints in local ignored folders use excluded alignment and BGC-level features and are not paper-model substitutes.

Treat this GitHub repository as the final analysis and manuscript release. The executable model artifact release is hosted separately on Hugging Face, and CHECKPOINTS.md records the source checkpoint hashes used to verify that release.

## Scientific Limits

- Labels are silver MIBiG-reference proxy groups, not direct chemical identity or pathway validation.
- The benchmark contains 16 eligible held-out families, limiting statistical power.
- Low Pfam overlap does not prove an alternative biochemical pathway or producer.
- The reported timings exclude Pfam and ESM preprocessing.
- The fixed split controls reference-group overlap, not all possible sequence homology.

## Citation and License

Citation metadata will be added after the preprint receives its persistent identifier. No reuse license has yet been assigned; the source is public for inspection, but reuse requires permission from the authors.
