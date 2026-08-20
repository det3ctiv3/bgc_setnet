"""Regenerate manuscript figures from the canonical final result tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LABELS = {
    "raw_esm_mean": "Raw ESM mean",
    "setnet": "Pfam-augmented SetNet",
    "pfam": "Pfam Jaccard",
    "ensemble": "ESM/Pfam-SetNet ensemble",
    "weighted_pfam": "Weighted Pfam Jaccard",
}
METRIC_LABELS = {
    "recall@50": "Recall@50",
    "mrr": "MRR",
    "map": "MAP",
    "ndcg@50": "nDCG@50",
}
ORDER = list(LABELS)
COLORS = ["#4C78A8", "#72B7B2", "#F2CF5B", "#E45756", "#B279A2"]


def load(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(root / "summary.csv")
    groups = pd.read_csv(root / "family_means_across_seeds_long.csv")
    required = {"method", "metric", "mean", "ci_lower", "ci_upper", "groups"}
    if not required.issubset(summary.columns):
        raise ValueError(f"Missing summary columns: {required - set(summary.columns)}")
    if not {"method", "group_id", "recall@50"}.issubset(groups.columns):
        raise ValueError("Family table is missing required columns")
    return summary, groups


def summary_figure(summary: pd.DataFrame, output: Path) -> None:
    frame = summary[summary.metric.eq("recall@50")].set_index("method").loc[ORDER].reset_index()
    x = np.arange(len(frame))
    yerr = np.vstack((frame["mean"] - frame["ci_lower"], frame["ci_upper"] - frame["mean"]))
    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    ax.bar(x, frame["mean"], yerr=yerr, capsize=4, color=COLORS, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels([LABELS[m] for m in frame.method], rotation=18, ha="right")
    ax.set_ylabel("Recall@50")
    ax.set_ylim(0.70, 0.96)
    ax.set_title("Group-disjoint held-out retrieval")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def metrics_figure(summary: pd.DataFrame, output: Path) -> None:
    metrics = ["recall@50", "mrr", "map", "ndcg@50"]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6), sharex=False)
    for ax, metric in zip(axes.flat, metrics):
        frame = summary[summary.metric.eq(metric)].set_index("method").loc[ORDER].reset_index()
        x = np.arange(len(frame))
        yerr = np.vstack((frame["mean"] - frame["ci_lower"], frame["ci_upper"] - frame["mean"]))
        ax.bar(x, frame["mean"], yerr=yerr, capsize=2, color=COLORS, edgecolor="black", linewidth=0.3)
        ax.set_title(METRIC_LABELS[metric])
        ax.set_xticks(x); ax.set_xticklabels([LABELS[m] for m in frame.method], rotation=18, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def family_scatter(groups: pd.DataFrame, output: Path) -> None:
    pivot = groups.pivot(index="group_id", columns="method", values="recall@50")
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), sharex=True, sharey=True)
    for ax, method, color in zip(axes, ["ensemble", "weighted_pfam"], ["#E45756", "#B279A2"]):
        frame = pivot.dropna(subset=["pfam", method])
        ax.scatter(frame["pfam"], frame[method], s=34, alpha=0.85, color=color, edgecolor="white", linewidth=0.4)
        bounds = [min(frame["pfam"].min(), frame[method].min()), max(frame["pfam"].max(), frame[method].max())]
        ax.plot(bounds, bounds, "--", color="0.35", linewidth=1)
        ax.set_title(LABELS[method])
        ax.set_xlabel("Pfam Jaccard Recall@50")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Method Recall@50")
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def distribution_figure(groups: pd.DataFrame, output: Path) -> None:
    values = [groups.loc[groups.method.eq(method), "recall@50"].to_numpy() for method in ORDER]
    fig, ax = plt.subplots(figsize=(8.2, 4.3))
    parts = ax.violinplot(values, showmeans=False, showmedians=True, showextrema=False)
    for body, color in zip(parts["bodies"], COLORS):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_alpha(0.75)
    ax.set_xticks(range(1, len(ORDER) + 1)); ax.set_xticklabels([LABELS[m] for m in ORDER], rotation=18, ha="right")
    ax.set_ylabel("Per-family Recall@50")
    ax.set_title("Distribution across eligible held-out families")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("results/dgx_final"))
    parser.add_argument("--output", type=Path, default=Path("paper"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summary, groups = load(args.data)
    summary_figure(summary, args.output / "fig_final_summary.pdf")
    metrics_figure(summary, args.output / "fig_final_metrics.pdf")
    family_scatter(groups, args.output / "fig_final_family_scatter.pdf")
    distribution_figure(groups, args.output / "fig_final_distribution.pdf")


if __name__ == "__main__":
    main()
