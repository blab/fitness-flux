#!/usr/bin/env python3
"""Build the similarity-accuracy component's data.json.

Two panels' worth of payload, from forecast_similarity.py's summaries:

  * the same forecasts scored by HA1 distance, in amino acids, against lead time;
  * the same forecasts scored by clade label, as mean absolute frequency error
    over the union of the two clade sets, against lead time.

The two panels carry different units and therefore separate axes; they are small
multiples, never a dual axis.

Takes one --track HORIZON SUMMARY DETAIL group per forecast horizon; the longest
horizon supplies the panels.
"""

import argparse
import json
from collections import defaultdict

import viz_io

# MLR blue / naive red, the paper's cool-warm mapping, validated against the
# light chart surface (dataviz skill: all six checks pass for this pair).
MODELS = [
    {"key": "mlr", "label": "MLR", "color": "#2a78d6"},
    {"key": "naive", "label": "Naïve", "color": "#e34948"},
]


def pct(value):
    return None if value is None else round(value * 100, 4)


def curve_for(summary):
    """Lead-time curves for both scores, from the summary's binned aggregate.

    HA1 distance stays in amino acids; the clade-label error is converted to
    percent to match the units the rest of the paper reports error in.
    """
    aa, mae = [], []
    for entry in summary["curve"]:
        lead = entry["lead"]
        aa.append({"lead": lead, "mlr": entry.get("mlr_eps0"), "naive": entry.get("naive_eps0")})
        mae.append({"lead": lead, "mlr": pct(entry.get("mlr_mae")), "naive": pct(entry.get("naive_mae"))})
    return aa, mae


def pairs_for(detail_path, bin_days, value_column, scale=1.0):
    """Faint per-window-pair curves, re-binned on the same weekly lead grid."""
    order = []
    bucket = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in viz_io.read_tsv(detail_path):
        label = f"{row['window']}->{row['next_window']}"
        if label not in bucket:
            order.append(label)
        lead = int(row["lead_days"])
        bucket[label][lead // bin_days][row["model"]].append(float(row[value_column]))

    out = []
    for label in order:
        points = []
        for b in sorted(bucket[label]):
            entry = {"lead": b * bin_days + (bin_days - 1) / 2.0}
            for model in ("mlr", "naive"):
                values = bucket[label][b].get(model, [])
                entry[model] = round(sum(values) / len(values) * scale, 4) if values else None
            points.append(entry)
        out.append({"label": label, "points": points})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--track",
        nargs=3,
        action="append",
        required=True,
        metavar=("HORIZON", "SUMMARY", "DETAIL"),
        help="one group per forecast horizon; repeat. The longest supplies the panels.",
    )
    parser.add_argument("--label", default="H3N2")
    parser.add_argument(
        "--gene-label",
        default="HA1",
        help="scored region, for the axis/metric labels (H3N2: HA1; SARS-CoV-2: spike S1)",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--meta-output",
        default=None,
        help="optional; the similarity components' meta.json is normally written "
        "once by the viz_similarity_meta rule",
    )
    args = parser.parse_args()

    tracks = []
    for horizon, summary_path, detail_path in args.track:
        with open(summary_path) as handle:
            summary = json.load(handle)
        tracks.append({
            "horizon": int(horizon),
            "summary": summary,
            "detail": detail_path,
        })
    tracks.sort(key=lambda t: t["horizon"])

    primary = max(tracks, key=lambda t: t["horizon"])
    summary = primary["summary"]
    aa_curve, mae_curve = curve_for(summary)
    bin_days = summary["lead_bin_days"]

    forecast = summary["overall"]["forecast"]
    metrics = [
        {
            "key": "aa",
            "label": f"Scored by {args.gene_label} distance",
            "ylabel": f"Forecast error ({args.gene_label} amino acids)",
            "curve": aa_curve,
            "pairs": pairs_for(primary["detail"], bin_days, "error_aa"),
            "mlr_mean": round(forecast["mlr_eps0"], 4),
            "naive_mean": round(forecast["naive_eps0"], 4),
            "unit": "aa",
        },
        {
            "key": "mae",
            "label": "Scored by clade label",
            "ylabel": "Mean absolute frequency error (%)",
            "curve": mae_curve,
            "pairs": pairs_for(primary["detail"], bin_days, "mae", scale=100.0),
            "mlr_mean": pct(forecast["mlr_mae"]),
            "naive_mean": pct(forecast["naive_mae"]),
            "unit": "%",
        },
    ]
    for metric in metrics:
        naive, mlr = metric["naive_mean"], metric["mlr_mean"]
        metric["gap_pct"] = round((naive - mlr) / naive * 100, 2) if naive else None

    viz_io.write_json(
        args.output,
        {
            "label": args.label,
            "horizon": primary["horizon"],
            "n_pairs": summary["n_pairs"],
            "max_distance": summary["max_distance"],
            "n_clades": summary["n_clades"],
            "n_profile_columns": summary["n_profile_columns"],
            "metrics": metrics,
            "endpoints": summary.get("endpoints", []),
            "models": MODELS,
        },
    )
    if args.meta_output:
        viz_io.write_json(
            args.meta_output,
            {"datasets": [{"id": "h3n2_clades", "label": "Similarity-aware accuracy"}], "default": "h3n2_clades"},
            indent=2,
        )


if __name__ == "__main__":
    main()
