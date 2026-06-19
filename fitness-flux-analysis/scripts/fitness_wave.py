#!/usr/bin/env python3
"""Characterize the fitness wave (the "flux") through time.

Port of the "Characterizing fitness wave" section of ``fitness-flux.nb``. Using
the scaffolded log fitness per variant and the per-date normalized frequencies:

  * location(t)  = sum_v freq_v(t) * logfit_v        (population mean log fitness,
                                                       the centre of the wave)
  * variance(t)  = sum_v freq_v(t) * (logfit_v - location(t))^2
  * velocity(t)  = generationTime * dLocation/dt      (per-generation fitness flux),
                   a centred finite difference over a fixed-day window.

Fisher's fundamental theorem predicts velocity ~ variance; the variance-vs-velocity
regression slope (near 1) is the scientific punchline, written to the summary JSON.
"""

import argparse
import csv
import json
import math

import numpy as np

import ff_io

DAYS_PER_YEAR = 365.0


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="e.g. sarscov2_clades")
    parser.add_argument("--scaffolded", required=True)
    parser.add_argument("--frequencies", required=True)
    parser.add_argument(
        "--generation-time",
        type=float,
        default=3.2,
        help="generation time in days (velocity is scaled to per-generation)",
    )
    parser.add_argument(
        "--velocity-window",
        type=int,
        default=None,
        help="finite-difference window in days (default 14 for sarscov2, else 60)",
    )
    parser.add_argument("--timeseries-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def read_scaffolded(path):
    with open(path) as handle:
        return {
            row["variant"]: float(row["log_fitness"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_frequencies(path):
    """Return ordered dates and {date: {variant: frequency}}."""
    by_date = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            date = row["date"]
            by_date.setdefault(date, {})[row["variant"]] = float(row["frequency"])
    return sorted(by_date), by_date


def location_and_variance(dates, by_date, logfit):
    location = []
    variance = []
    for date in dates:
        freqs = by_date[date]
        loc = sum(freq * logfit[v] for v, freq in freqs.items() if v in logfit)
        var = sum(
            freq * (logfit[v] - loc) ** 2 for v, freq in freqs.items() if v in logfit
        )
        location.append(loc)
        variance.append(var)
    return location, variance


def velocity_series(dates, location, generation_time, window):
    gen_years = generation_time / DAYS_PER_YEAR
    series = []  # (midpoint_date, velocity)
    decimals = [ff_io.decimal_year(d) for d in dates]
    for i in range(window, len(dates)):
        years = decimals[i] - decimals[i - window]
        if years <= 0:
            continue
        distance = location[i] - location[i - window]
        velocity = gen_years * distance / years
        series.append((dates[i - window // 2], velocity))
    return series


def main():
    args = parse_args()
    window = args.velocity_window
    if window is None:
        window = 14 if args.dataset.startswith("sarscov2") else 60

    logfit = read_scaffolded(args.scaffolded)
    dates, by_date = read_frequencies(args.frequencies)
    location, variance = location_and_variance(dates, by_date, logfit)
    velocity = velocity_series(dates, location, args.generation_time, window)
    velocity_by_date = dict(velocity)

    with open(args.timeseries_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["date", "decimal_date", "mean_log_fitness", "variance", "velocity"]
        )
        for date, loc, var in zip(dates, location, variance):
            vel = velocity_by_date.get(date)
            writer.writerow(
                [
                    date,
                    f"{ff_io.decimal_year(date):.4f}",
                    f"{loc:.6f}",
                    f"{var:.8f}",
                    "" if vel is None else f"{vel:.8f}",
                ]
            )
    ff_io.log(f"Wrote {len(dates)} flux-timeseries rows to {args.timeseries_output}")

    # Fisher check: regress velocity (at midpoint dates) on variance at those dates
    var_by_date = dict(zip(dates, variance))
    paired = [
        (var_by_date[d], v) for d, v in velocity if d in var_by_date
    ]
    summary = {
        "dataset": args.dataset,
        "generation_time_days": args.generation_time,
        "velocity_window_days": window,
        "avg_variance": float(np.mean(variance)),
        "avg_sd": float(np.mean([math.sqrt(v) for v in variance])),
        "avg_velocity": float(np.mean([v for _, v in velocity])) if velocity else None,
    }
    decimals = [ff_io.decimal_year(d) for d in dates]
    total_years = decimals[-1] - decimals[0]
    if total_years > 0:
        summary["simple_velocity"] = (
            (args.generation_time / DAYS_PER_YEAR)
            * (location[-1] - location[0])
            / total_years
        )
    if len(paired) >= 2:
        xs = np.array([p[0] for p in paired])
        ys = np.array([p[1] for p in paired])
        slope, intercept = np.polyfit(xs, ys, 1)
        r = float(np.corrcoef(xs, ys)[0, 1])
        summary["variance_vs_velocity"] = {
            "slope": float(slope),
            "intercept": float(intercept),
            "r_squared": r * r,
            "n": len(paired),
        }
    with open(args.summary_output, "w") as handle:
        json.dump(summary, handle, indent=2)
    ff_io.log(f"Wrote flux summary to {args.summary_output}")


if __name__ == "__main__":
    main()
