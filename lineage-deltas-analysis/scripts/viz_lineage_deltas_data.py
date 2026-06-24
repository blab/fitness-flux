#!/usr/bin/env python3
"""Merge the branch / predictor / ESM delta tables into one JSON for the
``lineage-deltas`` press component.

Each parent->child lineage branch (per season) becomes one point carrying its
change in log fitness plus every per-branch "predictor" delta (mutation-region
counts, mutation-effect scores, ESM embedding distances) under a single
friendly-keyed ``values`` map. The component picks which predictors to show as
panels via the ``predictors=`` figure attribute and subsets points by ``period``
(early/late) for its time toggle, so all the joining lives here.
"""

import argparse
import csv
import json

# Friendly predictor key -> (source column in branch_deltas.tsv, display label).
REGION_KEYS = {
    "nuc": ("delta_nuc_muts", "All nucleotide"),
    "spike": ("delta_spike_muts", "Spike"),
    "s1": ("delta_s1_muts", "Spike S1"),
    "s2": ("delta_s2_muts", "Spike S2"),
    "ntd": ("delta_ntd_muts", "Spike NTD"),
    "rbd": ("delta_rbd_muts", "Spike RBD"),
    "orf1ab": ("delta_orf1ab_muts", "ORF1ab"),
    "accessory": ("delta_accessory_muts", "Accessory"),
    "nonspike": ("delta_nonspike_muts", "Non-spike"),
}

# Friendly key (== predictor name in predictor_deltas.tsv) -> display label.
PREDICTOR_LABELS = {
    "covfit": "CoVFit",
    "evescape": "EvEscape",
    "dms_ace2_pseudovirus": "DMS ACE2 (pseudovirus)",
    "dms_entry_pseudovirus": "DMS entry (pseudovirus)",
    "dms_escape_pseudovirus": "DMS escape (pseudovirus)",
    "dms_ace2_yeast": "DMS ACE2 (yeast)",
    "dms_escape_yeast": "DMS escape (yeast)",
    "dms_expression_yeast": "DMS expression (yeast)",
    "dms_evescape": "DMS EVEscape",
    "dms_spike_muts": "DMS spike muts",
}

# ESM embedding columns in esm_deltas.tsv -> display label.
ESM_LABELS = {
    "esm_650M_pretrained": "ESM 650M pretrained",
    "esm_650M_fine_tuned": "ESM 650M fine-tuned",
    "esm_3B_pretrained": "ESM 3B pretrained",
    "esm_3B_fine_tuned": "ESM 3B fine-tuned",
}

# Early = Jan 2020 - Jun 2022 windows; everything else is late.
EARLY_SEASONS = ["2020", "2020-21", "2021", "2021-22"]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--branch-deltas", required=True)
    parser.add_argument("--predictor-deltas", required=True)
    parser.add_argument("--esm-deltas", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def branch_key(row):
    return (row["timepoint"], row["parent"], row["child"])


def read_branches(path):
    """(timepoint, parent, child) -> {timepoint, parent, child, delta_log_fitness, values{region keys}}."""
    points = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values = {}
            for key, (column, _label) in REGION_KEYS.items():
                v = to_float(row.get(column))
                if v is not None:
                    values[key] = v
            points[branch_key(row)] = {
                "timepoint": row["timepoint"],
                "parent": row["parent"],
                "child": row["child"],
                "delta_log_fitness": to_float(row["delta_log_fitness"]),
                "values": values,
            }
    return points


def attach_predictors(points, path):
    """Long predictor_deltas.tsv -> merge delta_predictor onto matching branches."""
    seen = set()
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            point = points.get(branch_key(row))
            if point is None:
                continue
            v = to_float(row["delta_predictor"])
            if v is not None:
                point["values"][row["predictor"]] = v
                seen.add(row["predictor"])
    return seen


def attach_esm(points, path):
    """Wide esm_deltas.tsv -> merge each ESM column onto matching branches."""
    seen = set()
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            point = points.get(branch_key(row))
            if point is None:
                continue
            for column in ESM_LABELS:
                v = to_float(row.get(column))
                if v is not None:
                    point["values"][column] = v
                    seen.add(column)
    return seen


def build_predictor_manifest(seen_predictors, seen_esm):
    manifest = {}
    for key, (_column, label) in REGION_KEYS.items():
        manifest[key] = {"label": label, "group": "Mutation region"}
    for key in seen_predictors:
        manifest[key] = {
            "label": PREDICTOR_LABELS.get(key, key),
            "group": "Effect predictor",
        }
    for key in seen_esm:
        manifest[key] = {"label": ESM_LABELS[key], "group": "ESM embedding"}
    return manifest


def main():
    args = parse_args()
    points = read_branches(args.branch_deltas)
    seen_predictors = attach_predictors(points, args.predictor_deltas)
    seen_esm = attach_esm(points, args.esm_deltas)

    ordered_points = sorted(
        points.values(), key=lambda p: (p["timepoint"], p["parent"], p["child"])
    )
    seasons = sorted({p["timepoint"] for p in ordered_points})
    early = [s for s in seasons if s in EARLY_SEASONS]
    for point in ordered_points:
        point["period"] = "early" if point["timepoint"] in EARLY_SEASONS else "late"

    document = {
        "dataset": "sarscov2_lineages",
        "seasons": seasons,
        "early": early,
        "predictors": build_predictor_manifest(seen_predictors, seen_esm),
        "points": ordered_points,
    }

    with open(args.output, "w") as handle:
        json.dump(document, handle, indent=0)
        handle.write("\n")

    n_early = sum(1 for p in ordered_points if p["period"] == "early")
    print(
        f"Wrote {len(ordered_points)} points "
        f"({n_early} early / {len(ordered_points) - n_early} late) "
        f"across {len(seasons)} seasons and {len(document['predictors'])} predictors "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
