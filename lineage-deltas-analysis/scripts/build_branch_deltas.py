#!/usr/bin/env python3
"""Build per-branch changes in mutation counts and fitness.

Port of the "Branches" section of ``lineage-deltas.nb``. For each parent->child
lineage edge present in a season (with MLR fitness on both endpoints and
mutation counts on both endpoints), record the change in mutation count per gene
region and the change in per-season log fitness. The same edge may recur across
seasons, contributing one row per season it appears in.

This is the core unit of the analysis: the central mutation-vs-fitness scatter,
the spike / non-spike split, the slope-through-time and the linear model are all
derived from these branch deltas.
"""

import argparse
import csv

import ld_io

GENE_COLUMNS = [
    "nuc_muts",
    "spike_muts",
    "s1_muts",
    "s2_muts",
    "ntd_muts",
    "rbd_muts",
    "orf1ab_muts",
    "accessory_muts",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument("--seqcounts-dir", default="sequence-counts")
    parser.add_argument("--mut-counts", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_mut_counts(path):
    counts = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lineage = row["lineage"]
            counts[lineage] = {}
            for column in GENE_COLUMNS:
                try:
                    counts[lineage][column] = float(row[column])
                except (ValueError, KeyError):
                    counts[lineage][column] = None
    return counts


def main():
    args = parse_args()
    mut_counts = read_mut_counts(args.mut_counts)
    branches = ld_io.branches(args.mlr_dir, args.seqcounts_dir)

    out_columns = (
        [f"delta_{name}" for name in GENE_COLUMNS]
        + ["delta_nonspike_muts", "delta_s1_nonrbd_muts"]
    )
    rows = []
    for branch in branches:
        child, parent = branch["child"], branch["parent"]
        if child not in mut_counts or parent not in mut_counts:
            continue
        deltas = {}
        ok = True
        for column in GENE_COLUMNS:
            cv, pv = mut_counts[child][column], mut_counts[parent][column]
            if cv is None or pv is None:
                ok = False
                break
            deltas[f"delta_{column}"] = cv - pv
        if not ok:
            continue
        deltas["delta_nonspike_muts"] = (
            deltas["delta_nuc_muts"] - deltas["delta_spike_muts"]
        )
        # RBD (319-541) is nested within S1 (14-685), so S1-outside-RBD is exact.
        deltas["delta_s1_nonrbd_muts"] = (
            deltas["delta_s1_muts"] - deltas["delta_rbd_muts"]
        )
        rows.append({**branch, **deltas})

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        header = ["timepoint", "parent", "child", "delta_log_fitness"] + out_columns
        writer.writerow(header)
        for row in rows:
            writer.writerow(
                [row["timepoint"], row["parent"], row["child"]]
                + [f"{row['delta_log_fitness']:.6f}"]
                + [f"{row[col]:g}" for col in out_columns]
            )
    ld_io.log(f"Wrote {len(rows)} branch-delta rows to {args.output}")


if __name__ == "__main__":
    main()
