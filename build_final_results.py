"""Build canonical manuscript tables from the final leakage-free DGX runs."""

from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


METRICS = ["recall@50", "mrr", "map", "ndcg@50"]
SEED_RE = re.compile(r"seed-(\d+)-evaluation$")


def exact_signed_rank(delta: np.ndarray) -> tuple[float, float, int]:
    delta = np.asarray(delta, dtype=float)
    delta = delta[~np.isclose(delta, 0.0, atol=1e-12)]
    ranks = rankdata(np.abs(delta), method="average")
    observed = min(float(ranks[delta > 0].sum()), float(ranks[delta < 0].sum()))
    count = 0
    for signs in itertools.product((-1, 1), repeat=len(delta)):
        signs = np.asarray(signs)
        statistic = min(float(ranks[signs > 0].sum()), float(ranks[signs < 0].sum()))
        count += statistic <= observed + 1e-12
    return observed, count / (2 ** len(delta)), len(delta)


def bootstrap_interval(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    samples = rng.choice(values, size=(20_000, len(values)), replace=True).mean(axis=1)
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def load_method(root: Path, roles: dict[str, str], selected: bool = False) -> pd.DataFrame:
    records = []
    for run_dir in sorted(root.glob("paper-*-evaluation")):
        match = SEED_RE.search(run_dir.name)
        if not match:
            continue
        groups = pd.read_csv(run_dir / "group_results.csv")
        metadata = json.loads((run_dir / "metadata.json").read_text())
        if selected:
            method = f"ensemble_validation_alpha_{float(metadata['selected_alpha']):.1f}"
            groups = groups[groups["method"].eq(method)].copy()
            groups["method"] = "ensemble"
        else:
            groups = groups[groups["method"].isin(roles)].copy()
            groups["method"] = groups["method"].map(roles)
        groups["seed"] = int(match.group(1))
        records.append(groups[["method", "group_id", "seed", *METRICS, "tie_fraction"]])
    if not records:
        raise ValueError(f"No result files found under {root}")
    return pd.concat(records, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pfam-setnet", type=Path, default=Path("results/dgx_pfam_setnet/evaluations"))
    parser.add_argument("--weighted-pfam", type=Path, default=Path("results/dgx_weighted_pfam/evaluations"))
    parser.add_argument("--output", type=Path, default=Path("results/dgx_final"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    pfam_setnet = load_method(
        args.pfam_setnet,
        {"raw_esm_mean": "raw_esm_mean", "setnet": "setnet", "pfam_jaccard_max": "pfam"},
    )
    pfam_setnet = pd.concat(
        [pfam_setnet, load_method(args.pfam_setnet, {}, selected=True)], ignore_index=True
    )
    weighted = load_method(
        args.weighted_pfam,
        {"pfam_jaccard_max": "pfam", "weighted_pfam_jaccard": "weighted_pfam"},
    )
    weighted = weighted[weighted["method"].eq("weighted_pfam")]
    all_results = pd.concat([pfam_setnet, weighted], ignore_index=True)
    all_results.to_csv(args.output / "per_seed_group_results.csv", index=False)

    family = all_results.groupby(["method", "group_id"])[METRICS + ["tie_fraction"]].mean().reset_index()
    family.to_csv(args.output / "family_means_across_seeds_long.csv", index=False)
    wide = family.pivot(index="group_id", columns="method", values=METRICS)
    wide.columns = [f"{method}__{metric}" for metric, method in wide.columns]
    wide.reset_index().to_csv(args.output / "family_means_across_seeds.csv", index=False)

    rng = np.random.default_rng(20260815)
    summary_rows = []
    for method in ["raw_esm_mean", "setnet", "pfam", "ensemble", "weighted_pfam"]:
        frame = family[family["method"].eq(method)]
        for metric in METRICS:
            values = frame[metric].to_numpy(dtype=float)
            low, high = bootstrap_interval(values, rng)
            summary_rows.append({
                "method": method, "metric": metric, "mean": float(values.mean()),
                "ci_lower": low, "ci_upper": high, "groups": len(values),
            })
    pd.DataFrame(summary_rows).to_csv(args.output / "summary.csv", index=False)

    comparisons = []
    for method in ["ensemble", "weighted_pfam"]:
        for metric in METRICS:
            left = family[family["method"].eq(method)].set_index("group_id")[metric]
            right = family[family["method"].eq("pfam")].set_index("group_id")[metric]
            common = left.index.intersection(right.index)
            delta = (left.loc[common] - right.loc[common]).to_numpy()
            statistic, p_value, nonzero = exact_signed_rank(delta)
            comparisons.append({
                "comparison": f"{method}_vs_pfam", "method": method, "metric": metric,
                "n_families": len(delta), "n_nonzero": nonzero,
                "mean_delta": float(delta.mean()), "median_delta": float(np.median(delta)),
                "wins": int((delta > 0).sum()), "ties": int((delta == 0).sum()),
                "losses": int((delta < 0).sum()), "wilcoxon_statistic": statistic,
                "p_value": p_value,
            })
    comparisons = pd.DataFrame(comparisons)
    primary = comparisons[comparisons["metric"].eq("recall@50")].sort_values("p_value")
    running = 0.0
    holm = {}
    for rank, (index, row) in enumerate(primary.iterrows()):
        running = max(running, min(1.0, float(row["p_value"]) * (len(primary) - rank)))
        holm[index] = running
    comparisons["holm_p_value"] = [holm.get(index, np.nan) for index in comparisons.index]
    comparisons.to_csv(args.output / "paired_comparisons.csv", index=False)
    json.dump({"seeds": [20260810, 20260811, 20260812, 20260813, 20260814],
               "eligible_test_groups": 16, "split_sha256": "dc26fae17e54fd2ad41a9e10353b3da3e0aacf3144b64f6ee62e8341b4360555",
               "unit_of_inference": "family mean across five seeds", "unknown_candidates": True},
              (args.output / "metadata.json").open("w"), indent=2)
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
