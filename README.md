# BGC-SetNet Retrieval Benchmark

Analysis code, result artifacts, and manuscript source for:

> Sequence-Derived Representations versus Pfam-Domain Content for Biosynthetic Gene Cluster Retrieval

This repository reports a group-disjoint, five-seed evaluation of biosynthetic gene cluster (BGC) retrieval in a *Streptomyces griseus* species-family atlas. The central result is negative but informative: explicit Pfam-domain content remains the strongest signal on this benchmark.

## Final Results

The effective inference unit is the held-out family (n = 16). Values below are five-seed family means.

| Method | Recall@50 | MRR | MAP | nDCG@50 |
| --- | ---: | ---: | ---: | ---: |
| Raw ESM mean | 0.7946 | 0.2550 | 0.7251 | 0.8078 |
| Pfam-augmented BGC-SetNet | 0.8472 | 0.2786 | 0.7771 | 0.8502 |
| Pfam Jaccard | 0.8788 | 0.3071 | 0.8480 | 0.9042 |
| ESM/Pfam-SetNet ensemble | 0.8769 | 0.3096 | 0.8503 | 0.9058 |
| Weighted Pfam Jaccard | 0.8789 | 0.3069 | 0.8477 | 0.9040 |

Weighted Pfam exceeds unweighted Pfam by only 0.00003 Recall@50 (p = 0.50). The ESM/Pfam-SetNet ensemble is lower by 0.00198 (p = 0.625). Neither difference supports an improvement claim.

## Release Scope

This release contains:

- build_final_results.py: aggregation, bootstrap intervals, and exact paired tests.
- plot_final_results.py: deterministic regeneration of the four manuscript figures.
- results/dgx_final/: compact canonical aggregate and per-seed family results.
- paper/main.tex: canonical manuscript source.
- paper/fig_final_*.pdf: figures generated from the canonical result tables.
- CHECKPOINTS.md: expected checkpoint hashes and model-release status.

The original data-preparation scripts remain in the repository. The obsolete July model checkpoints and local large assets are intentionally excluded.

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

Build the manuscript from the canonical paper directory:

~~~bash
cd paper
tectonic main.tex
~~~

## Reproducibility Boundary

The exact v5 Pfam-augmented SetNet checkpoints, v6 weighted-Pfam checkpoints, weighted-Pfam training implementation, training-only Pfam vocabulary, and frozen split manifest are not present on this workstation. They must be recovered from the original DGX run archive before a paper-model release can be made on Hugging Face. The available July checkpoints use excluded alignment and BGC-level features and are not substitutes.

Until those files are recovered and hash-verified, treat this repository as the final analysis and manuscript release, not as a complete executable training package. See CHECKPOINTS.md for the required SHA-256 values.

## Scientific Limits

- Labels are silver MIBiG-reference proxy groups, not direct chemical identity or pathway validation.
- The benchmark contains 16 eligible held-out families, limiting statistical power.
- Low Pfam overlap does not prove an alternative biochemical pathway or producer.
- The reported timings exclude Pfam and ESM preprocessing.
- The fixed split controls reference-group overlap, not all possible sequence homology.

## Citation and License

Citation metadata will be added after the preprint receives its persistent identifier. No reuse license has yet been assigned; the source is public for inspection, but reuse requires permission from the authors.
