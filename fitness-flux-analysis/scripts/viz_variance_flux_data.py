#!/usr/bin/env python3
"""
Build the variance-vs-flux component's data.json: the daily fitness-wave
timeseries (population fitness variance and flux/velocity per day) plus the
precomputed variance-vs-velocity regression. The component renders the four-panel
"fitness wave" figure (variance timeseries, variance-vs-velocity scatter,
velocity timeseries, per-year velocity means) from these daily points.

Emits per record {date, variance, velocity}; velocity is null for the early dates
that precede the velocity window (kept so the variance timeseries stays complete).
Mirrors viz_frequency_panels_data.py. Dates stay ISO strings; the component
constructs Date objects.
"""
import argparse
import json

import viz_io


def build_points(timeseries_path):
    points = []
    for row in viz_io.read_tsv(timeseries_path):
        variance = viz_io.num(row["variance"])
        if variance is None:
            continue
        points.append(
            {
                "date": row["date"],
                "variance": variance,
                "velocity": viz_io.num(row["velocity"]),
            }
        )
    return points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeseries", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.summary) as handle:
        summary = json.load(handle)

    viz_io.write_json(
        args.output,
        {
            "points": build_points(args.timeseries),
            "fit": summary.get("variance_vs_velocity", {}),
        },
    )


if __name__ == "__main__":
    main()
