"""
ESM-2 Inference: Embed all 294,604 proteins into 1280-d vectors.
Saves to HDF5 with gene_id as key for fast random access during training.

Hardware: NVIDIA DGX Spark (128GB unified memory)
Model: esm2_t33_650M_UR50D (650M params, 33 layers, 1280-d output)
Time estimate: ~6-10 hours at batch_size=16

Usage:
    python esm2_inference.py
    python esm2_inference.py --batch_size 8    # if OOM, reduce batch size
    python esm2_inference.py --resume          # resume from last checkpoint
"""
import os
import argparse
import time
import json
import torch
import h5py
import numpy as np
from pathlib import Path

FASTA_PATH = "data/filtered_proteins.fasta"
OUTPUT_PATH = "data/esm2_embeddings.h5"
CHECKPOINT_PATH = "data/esm2_progress.json"
MAX_SEQ_LEN = 1022  # ESM-2 max tokens (excluding BOS/EOS)
EMBED_DIM = 1280


def parse_fasta(fasta_path):
    """Parse FASTA file into list of (gene_id, sequence) tuples."""
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(current_seq)))
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences.append((current_id, "".join(current_seq)))

    return sequences


def batch_by_length(sequences, batch_size):
    """Sort by length and batch — minimizes padding waste."""
    sorted_seqs = sorted(sequences, key=lambda x: len(x[1]))
    batches = []
    for i in range(0, len(sorted_seqs), batch_size):
        batches.append(sorted_seqs[i:i + batch_size])
    return batches


def load_model():
    """Load ESM-2 650M model and alphabet."""
    print("Loading ESM-2 model (esm2_t33_650M_UR50D)...")
    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
    model = model.eval()

    if torch.cuda.is_available():
        model = model.cuda()
        print(f"  Running on GPU: {torch.cuda.get_device_name()}")
        print(f"  GPU memory: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")
    else:
        print("  WARNING: No GPU detected. This will be very slow.")

    model = model.half()  # fp16 for speed + memory
    print("  Using fp16 precision")

    return model, alphabet


def extract_embedding(model, alphabet, batch_sequences, device):
    """
    Run ESM-2 on a batch of sequences.
    Returns mean-pooled last-layer representations (1280-d per sequence).
    """
    batch_converter = alphabet.get_batch_converter()

    # Prepare data: list of (label, sequence) tuples
    data = [(gene_id, seq[:MAX_SEQ_LEN]) for gene_id, seq in batch_sequences]
    _, _, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad(), torch.cuda.amp.autocast():
        results = model(batch_tokens, repr_layers=[33], return_contacts=False)

    # Extract last layer representations
    token_reps = results["representations"][33]  # (batch, seq_len, 1280)

    # Mean pool over sequence (excluding BOS and EOS tokens)
    embeddings = []
    for i, (gene_id, seq) in enumerate(data):
        seq_len = min(len(seq), MAX_SEQ_LEN)
        # tokens: [BOS, aa1, aa2, ..., aaN, EOS, PAD, ...]
        # We want positions 1 to seq_len (inclusive) — the actual amino acids
        emb = token_reps[i, 1:seq_len + 1, :].mean(dim=0)  # (1280,)
        embeddings.append((gene_id, emb.cpu().float().numpy()))

    return embeddings


def save_checkpoint(processed_count, total, elapsed):
    """Save progress for resume capability."""
    with open(CHECKPOINT_PATH, "w") as f:
        json.dump({
            "processed": processed_count,
            "total": total,
            "elapsed_seconds": elapsed
        }, f)


def load_checkpoint():
    """Load progress from previous run."""
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH) as f:
            return json.load(f)
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch_size", type=int, default=16,
                        help="Proteins per batch. Reduce if OOM (try 8 or 4)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from last checkpoint")
    args = parser.parse_args()

    # Parse FASTA
    print(f"Parsing {FASTA_PATH}...")
    sequences = parse_fasta(FASTA_PATH)
    print(f"  Total proteins: {len(sequences):,}")

    seq_lengths = [len(s) for _, s in sequences]
    print(f"  Length range: {min(seq_lengths)} - {max(seq_lengths)} AA")
    print(f"  Median length: {sorted(seq_lengths)[len(seq_lengths)//2]} AA")
    truncated = sum(1 for l in seq_lengths if l > MAX_SEQ_LEN)
    if truncated:
        print(f"  Truncated to {MAX_SEQ_LEN}: {truncated:,} sequences")

    # Check resume state
    already_done = set()
    if args.resume and os.path.exists(OUTPUT_PATH):
        with h5py.File(OUTPUT_PATH, "r") as h5:
            already_done = set(h5.keys())
        print(f"\n  Resuming: {len(already_done):,} already embedded, {len(sequences) - len(already_done):,} remaining")
        sequences = [(gid, seq) for gid, seq in sequences if gid not in already_done]

    if not sequences:
        print("All proteins already embedded. Nothing to do.")
        return

    # Batch by length
    batches = batch_by_length(sequences, args.batch_size)
    print(f"  Batches: {len(batches):,} (batch_size={args.batch_size})")

    # Load model
    model, alphabet = load_model()
    device = next(model.parameters()).device

    # Open HDF5 for writing
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    mode = "a" if args.resume and os.path.exists(OUTPUT_PATH) else "w"
    h5 = h5py.File(OUTPUT_PATH, mode)

    # Run inference
    total_processed = len(already_done)
    total_proteins = total_processed + len(sequences)
    start_time = time.time()
    last_print = start_time

    print(f"\nStarting inference...")
    print(f"{'='*60}")

    for batch_idx, batch in enumerate(batches):
        try:
            embeddings = extract_embedding(model, alphabet, batch, device)

            for gene_id, emb in embeddings:
                h5.create_dataset(gene_id, data=emb.astype(np.float16))
                total_processed += 1

        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                torch.cuda.empty_cache()
                print(f"\n  OOM on batch {batch_idx} (max_len={max(len(s) for _, s in batch)})")
                print(f"  Splitting batch and retrying...")
                # Process one by one for this batch
                for single in batch:
                    try:
                        embeddings = extract_embedding(model, alphabet, [single], device)
                        for gene_id, emb in embeddings:
                            h5.create_dataset(gene_id, data=emb.astype(np.float16))
                            total_processed += 1
                    except RuntimeError:
                        torch.cuda.empty_cache()
                        print(f"  SKIPPED (too long): {single[0]} ({len(single[1])} AA)")
            else:
                raise

        # Progress reporting every 30 seconds
        now = time.time()
        if now - last_print > 30:
            elapsed = now - start_time
            rate = (total_processed - len(already_done)) / elapsed
            remaining = (total_proteins - total_processed) / max(rate, 0.01)
            pct = total_processed / total_proteins * 100
            print(f"  [{pct:5.1f}%] {total_processed:,}/{total_proteins:,} | "
                  f"{rate:.1f} proteins/sec | ETA: {remaining/3600:.1f}h")
            last_print = now
            h5.flush()
            save_checkpoint(total_processed, total_proteins, elapsed)

    h5.close()
    elapsed = time.time() - start_time

    # Final stats
    print(f"{'='*60}")
    print(f"\nDone!")
    print(f"  Proteins embedded: {total_processed:,}")
    print(f"  Time: {elapsed/3600:.1f} hours")
    print(f"  Output: {OUTPUT_PATH}")
    print(f"  Size: {os.path.getsize(OUTPUT_PATH) / 1e9:.2f} GB")

    # Cleanup checkpoint
    if os.path.exists(CHECKPOINT_PATH):
        os.remove(CHECKPOINT_PATH)

    # Verify
    with h5py.File(OUTPUT_PATH, "r") as h5:
        n_keys = len(h5.keys())
        sample_key = list(h5.keys())[0]
        sample_shape = h5[sample_key].shape
        print(f"\n  Verification:")
        print(f"    Entries: {n_keys:,}")
        print(f"    Shape per protein: {sample_shape}")
        print(f"    Dtype: {h5[sample_key].dtype}")


if __name__ == "__main__":
    main()
