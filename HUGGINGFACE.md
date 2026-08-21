# Hugging Face Model Release

Public model repository: https://huggingface.co/whiteh4t/bgc-setnet

This release contains the verified final-paper model artifacts recovered from the original journal rebuild run archive:

- Five Pfam-augmented BGC-SetNet `model.safetensors` checkpoints for seeds 20260810-20260814.
- Five weighted-Pfam Jaccard `model.safetensors` checkpoints for seeds 20260810-20260814.
- Per-seed configs, Pfam vocabularies, training histories, and evaluation metadata.
- Aggregate manuscript result tables under `results/dgx_final/`.
- A source snapshot under `source/` and a provenance manifest under `provenance/checkpoint_manifest.json`.

All ten source checkpoints were verified against the expected SHA-256 values before conversion to safetensors. The public package excludes the original `.pt` training checkpoints because their metadata contains workstation-specific absolute paths; the safetensors files strict-load against the included source.

The model repository is public. License terms remain marked as `other`/pending in the model card, so add a formal license before encouraging broad reuse.
