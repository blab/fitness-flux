#!/usr/bin/env python3
"""Compute per-branch changes in mutation-effect predictors vs fitness.

Port of the EvEscape / CoVFit / DMS predictor sections of ``lineage-deltas.nb``.
Each predictor is a per-lineage score; along a parent->child branch its change is
``predictor[child] - predictor[parent]``, paired with the per-season change in log
fitness. Output is long format so the dashboard can draw one correlation panel
per predictor.

Predictors:
  * covfit                   CoVFit predicted fitness (seq_name -> fitness_mean)
  * evescape                 mean EvEscape per mode_lineage (evescape.org)
  * dms_*                    Bloom-lab clade phenotypes (pseudovirus + yeast display)
"""

import argparse
import csv

import ld_io

# clade_phenotypes.csv column headers -> predictor name
DMS_COLUMNS = {
    "dms_spike_muts": "number spike muts from Wuhan-Hu-1",
    "dms_escape_pseudovirus": "spike pseudovirus DMS human sera escape relative to XBB.1.5",
    "dms_ace2_pseudovirus": "spike pseudovirus DMS ACE2 binding relative to XBB.1.5",
    "dms_entry_pseudovirus": "spike pseudovirus DMS spike mediated entry relative to XBB.1.5",
    "dms_ace2_yeast": "RBD yeast-display DMS ACE2 affinity relative to XBB.1.5",
    "dms_expression_yeast": "RBD yeast-display DMS RBD expression relative to XBB.1.5",
    "dms_escape_yeast": "RBD yeast-display DMS escape relative to XBB.1.5",
    "dms_evescape": "EVEscape relative to XBB.1.5",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument("--seqcounts-dir", default="sequence-counts")
    parser.add_argument("--covfit", required=True)
    parser.add_argument("--evescape", required=True)
    parser.add_argument("--dms", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_evescape(path):
    """mode_lineage -> mean evescape across sequences."""
    totals = {}
    counts = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            try:
                value = float(row["evescape"])
            except (ValueError, KeyError, TypeError):
                continue
            lineage = row["mode_lineage"]
            totals[lineage] = totals.get(lineage, 0.0) + value
            counts[lineage] = counts.get(lineage, 0) + 1
    return {lineage: totals[lineage] / counts[lineage] for lineage in totals}


def load_predictors(args):
    predictors = {}
    predictors["covfit"] = ld_io.read_tsv_map(args.covfit, "seq_name", "fitness_mean")
    predictors["evescape"] = read_evescape(args.evescape)
    for name, column in DMS_COLUMNS.items():
        predictors[name] = ld_io.read_tsv_map(
            args.dms, "clade", column, delimiter=","
        )
    return predictors


def main():
    args = parse_args()
    predictors = load_predictors(args)
    branches = ld_io.branches(args.mlr_dir, args.seqcounts_dir)

    rows = []
    for name, scores in predictors.items():
        for branch in branches:
            child, parent = branch["child"], branch["parent"]
            if child in scores and parent in scores:
                rows.append(
                    {
                        "timepoint": branch["timepoint"],
                        "parent": parent,
                        "child": child,
                        "predictor": name,
                        "delta_predictor": scores[child] - scores[parent],
                        "delta_log_fitness": branch["delta_log_fitness"],
                    }
                )

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "timepoint",
                "parent",
                "child",
                "predictor",
                "delta_predictor",
                "delta_log_fitness",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["timepoint"],
                    row["parent"],
                    row["child"],
                    row["predictor"],
                    f"{row['delta_predictor']:.6f}",
                    f"{row['delta_log_fitness']:.6f}",
                ]
            )
    by_pred = {}
    for row in rows:
        by_pred[row["predictor"]] = by_pred.get(row["predictor"], 0) + 1
    ld_io.log(f"Wrote {len(rows)} predictor-delta rows to {args.output}: {by_pred}")


if __name__ == "__main__":
    main()
