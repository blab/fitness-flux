#!/usr/bin/env python3
"""
Build the forecast-accuracy component's data.json: the aggregate MAE(lead) curve
(MLR vs naive) plus faint per-window-pair curves, with absolute errors converted
to percent for the "%" y-axis. The component renders MLR-vs-naive mean absolute
error against forecast lead time (Abousamra et al. Fig 2B adaptation).

The aggregate `curve` comes straight from the summary JSON; the per-pair `points`
are re-binned here from the detail TSV using the same weekly lead bins. Mirrors
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--detail", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.summary) as handle:
        summary = json.load(handle)

    viz_io.write_json(
        args.output,
        {
            "curve": build_curve(summary),
            "pairs": build_pairs(args.detail, summary["lead_bin_days"]),
            "reference": pct(summary["reference"]),
            "models": MODELS,
        },
    )


if __name__ == "__main__":
    main()
