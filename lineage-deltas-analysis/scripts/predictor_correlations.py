#!/usr/bin/env python3
"""Summarize how well each predictor's branch delta tracks fitness change.

Port of the per-predictor correlation reporting in ``lineage-deltas.nb``
(EvEscape / CoVFit / DMS / ESM correlation sections). For each predictor we
report Pearson r, Spearman rho, regression slope and R^2 of its per-branch delta
against the per-branch change in log fitness, drawing from the branch-delta,
predictor-delta and ESM-delta tables.
"""

import argparse
import csv

import numpy as np
from scipy import stats

import ld_io

# mutation-count predictors live as columns in the branch-deltas table
MUTATION_PREDICTORS = [
    "delta_nuc_muts",
    "delta_spike_muts",
    "delta_nonspike_muts",
    "delta_rbd_muts",
    "delta_ntd_muts",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--branch-deltas", required=True)
    parser.add_argument("--predictor-deltas", required=True)
    parser.add_argument("--esm-deltas", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def correlate(xs, ys):
    xs = np.array(xs)
    ys = np.array(ys)
    if len(xs) < 3 or np.ptp(xs) == 0:
        return None
    pearson = float(stats.pearsonr(xs, ys)[0])
    spearman = float(stats.spearmanr(xs, ys)[0])
    slope, intercept = np.polyfit(xs, ys, 1)
    return {
        "pearson_r": pearson,
        "spearman_rho": spearman,
        "slope": float(slope),
        "r_squared": pearson * pearson,
        "n": len(xs),
    }


def collect(records, x_key, y_key="delta_log_fitness"):
    xs, ys = [], []
    for record in records:
        try:
            xs.append(float(record[x_key]))
            ys.append(float(record[y_key]))
        except (ValueError, KeyError, TypeError):
            continue
    return xs, ys


def read_rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    args = parse_args()
    results = []  # (predictor, stats)

    branch_rows = read_rows(args.branch_deltas)
    for column in MUTATION_PREDICTORS:
        xs, ys = collect(branch_rows, column)
        stat = correlate(xs, ys)
        if stat:
            results.append((column.replace("delta_", "").replace("_muts", ""), stat))

    predictor_rows = read_rows(args.predictor_deltas)
    by_predictor = {}
    for row in predictor_rows:
        by_predictor.setdefault(row["predictor"], []).append(row)
    for predictor, rows in by_predictor.items():
        xs, ys = collect(rows, "delta_predictor")
        stat = correlate(xs, ys)
        if stat:
            results.append((predictor, stat))

    esm_rows = read_rows(args.esm_deltas)
    for column in ["esm_650M_pretrained", "esm_650M_fine_tuned", "esm_3B_pretrained", "esm_3B_fine_tuned"]:
        xs, ys = collect(esm_rows, column)
        stat = correlate(xs, ys)
        if stat:
            results.append((column, stat))

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["predictor", "pearson_r", "spearman_rho", "slope", "r_squared", "n"]
        )
        for predictor, stat in results:
            writer.writerow(
                [
                    predictor,
                    f"{stat['pearson_r']:.4f}",
                    f"{stat['spearman_rho']:.4f}",
                    f"{stat['slope']:.6f}",
                    f"{stat['r_squared']:.4f}",
                    stat["n"],
                ]
            )
    ld_io.log(f"Wrote {len(results)} predictor correlations to {args.output}")


if __name__ == "__main__":
    main()
