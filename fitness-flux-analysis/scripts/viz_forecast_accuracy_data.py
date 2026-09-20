#!/usr/bin/env python3
"""
Build the forecast-accuracy component's combined data.json: one panel per virus
lineage, each carrying the aggregate MAE(lead) curve (MLR vs naive) plus faint
per-window-pair curves, with absolute errors converted to percent for the "%"
y-axis. The component renders these as a 2x2 grid of MLR-vs-naive mean absolute
error against forecast lead time (Abousamra et al. Fig 2B adaptation).

Takes one --panel KEY LABEL SUMMARY DETAIL group per lineage; panel order in the
output follows the argument order. Each panel's aggregate `curve` comes from the
summary JSON and its per-pair `points` are re-binned from the detail TSV using
the same weekly lead bins. Also writes the component's meta.json. Mirrors
viz_variance_flux_data.py.
"""
import argparse
import json
from collections import defaultdict

import viz_io

# Model color map: MLR blue, naive red (the paper's cool/warm mapping), using
# the validated categorical blue/red steps (dataviz skill; CVD-safe on the light
# surface) rather than the paper's darker navy/maroon.
MODELS = [
    {"key": "mlr", "label": "MLR", "color": "#2a78d6"},
    {"key": "naive", "label": "Naïve", "color": "#e34948"},
]


def pct(value):
    return None if value is None else round(value * 100, 4)


def build_curve(summary):
    return [
        {"lead": e["lead"], "mlr": pct(e["mlr"]), "naive": pct(e["naive"])}
        for e in summary["curve"]
    ]


def build_pairs(detail_path, bin_days):
    """Per-pair MAE(lead) curves, re-binned by the same weekly lead bins."""
    order = []
    bucket = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # label -> bin -> model -> [ae]
    for row in viz_io.read_tsv(detail_path):
        label = f"{row['window']}->{row['next_window']}"
        if label not in bucket:
            order.append(label)
        lead = int(row["lead_days"])
        bucket[label][lead // bin_days][row["model"]].append(float(row["abs_error"]))

    pairs = []
    for label in order:
        by_bin = bucket[label]
        points = []
        for b in sorted(by_bin):
            entry = {"lead": b * bin_days + (bin_days - 1) / 2.0}
            for model in ("mlr", "naive"):
                values = by_bin[b].get(model, [])
                entry[model] = pct(sum(values) / len(values)) if values else None
            points.append(entry)
        pairs.append({"label": label, "points": points})
    return pairs


def build_panel(key, label, summary_path, detail_path):
    with open(summary_path) as handle:
        summary = json.load(handle)
    return {
        "key": key,
        "label": label,
        "n_pairs": summary["n_pairs"],
        "curve": build_curve(summary),
        "pairs": build_pairs(detail_path, summary["lead_bin_days"]),
    }, pct(summary["reference"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        nargs=4,
        action="append",
        required=True,
        metavar=("KEY", "LABEL", "SUMMARY", "DETAIL"),
        help="one group per lineage; repeat. Panel order follows argument order.",
    )
    parser.add_argument("--output", required=True, help="combined data/all.json")
    parser.add_argument("--meta-output", required=True, help="component meta.json")
    args = parser.parse_args()

    panels = []
    reference = 5.0
    for key, label, summary_path, detail_path in args.panel:
        panel, reference = build_panel(key, label, summary_path, detail_path)
        panels.append(panel)

    viz_io.write_json(
        args.output,
        {"panels": panels, "reference": reference, "models": MODELS},
    )
    viz_io.write_json(
        args.meta_output,
        {"datasets": [{"id": "all", "label": "Forecast accuracy"}], "default": "all"},
        indent=2,
    )


if __name__ == "__main__":
    main()
