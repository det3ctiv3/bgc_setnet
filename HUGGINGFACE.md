# Hugging Face Model Release Gate

The paper model must not be published from the checkpoint files currently on this workstation.

The only local weights are July legacy checkpoints. They load successfully, but they use 1,284-dimensional gene inputs plus 13 BGC features, including alignment and target-derived fields excluded by the final manuscript protocol. Their hashes do not match the reported v5/v6 runs.

## Target Repository

Planned namespace: whiteh4t/bgc-setnet

Do not create or advertise this repository as the paper model until every required artifact in CHECKPOINTS.md has been recovered and verified.

## Required Contents

- README.md model card with intended use, limitations, silver-label warning, and citation.
- config.json for the final 1,280-to-256 SetNet encoder.
- Exact v5/v6 training and inference source.
- Five Pfam-augmented SetNet checkpoints, one per reported seed.
- Five weighted-Pfam checkpoints, one per reported seed.
- Training-only 2,218-token Pfam vocabulary with PAD/UNK mapping.
- Per-seed validation-selected ensemble weights: 0.1, 0.1, 0.1, 0.2, 0.4.
- Preprocessing specification for ESM-2 model revision, layer, pooling, gene order, relative positions, and padding.
- Frozen split manifest and SHA-256 verifier.
- Checkpoint manifest containing seed, split hash, validation-selection metadata, code commit, and file hash.
- A smoke-tested inference example.
- Explicit code, model, and data licenses.

## Upload Validation

Before upload:

1. Verify all ten SHA-256 values in CHECKPOINTS.md.
2. Strict-load each checkpoint with the recovered final implementation.
3. Run a finite-output and unit-normalization inference smoke test for every SetNet seed.
4. Reproduce the committed family-level summary from the recovered artifacts.
5. Convert tensor weights to safetensors where possible without discarding required metadata.
6. Upload all five seeds; do not select a seed using held-out test performance.

Large atlas, FASTA, and ESM embedding files should be released separately as a dataset only after provenance and redistribution rights are documented.
