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
import os
import sys

import numpy as np

import ff_io

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from generation_time import generation_time_for, load_aliasor  # noqa: E402

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
        "--generation-time-pre-omicron",
        type=float,
        default=None,
        help="if set, the per-timepoint generation time is a frequency-weighted mean "
        + "of per-variant tau (pre-Omicron variants use this, others --generation-time)",
    )
    parser.add_argument(
        "--variant-classification",
        choices=["clades", "lineages"],
        default=None,
        help="how to read variant names for the pre/post-Omicron split",
    )
    parser.add_argument(
        "--aliasing",
        type=str,
        default=None,
        help="optional local Pango alias_key.json for lineage classification",
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


def per_variant_tau(variants, tau_post, tau_pre, variant_classification, aliasor):
    """variant -> generation time. Uniform tau_post when no pre-Omicron split."""
    if tau_pre is None or variant_classification is None:
        return {v: tau_post for v in variants}
    return {
        v: generation_time_for(v, variant_classification, tau_pre, tau_post, aliasor)
        for v in variants
    }


def tau_bar_by_date(dates, by_date, tau_v, tau_post):
    """Frequency-weighted mean generation time at each date, tau_bar(t) = sum_v
    freq_v(t) * tau_v; the effective generation clock as the population shifts from
    pre-Omicron (5.0) to Omicron (3.2)."""
    out = {}
    for date in dates:
        freqs = by_date[date]
        total = sum(freqs.values())
        if total <= 0:
            out[date] = tau_post
        else:
            out[date] = sum(f * tau_v.get(v, tau_post) for v, f in freqs.items()) / total
    return out


def velocity_series(dates, location, tau_bar, window):
    series = []  # (midpoint_date, velocity)
    decimals = [ff_io.decimal_year(d) for d in dates]
    for i in range(window, len(dates)):
        years = decimals[i] - decimals[i - window]
        if years <= 0:
            continue
        distance = location[i] - location[i - window]
        midpoint = dates[i - window // 2]
        gen_years = tau_bar[midpoint] / DAYS_PER_YEAR
        velocity = gen_years * distance / years
        series.append((midpoint, velocity))
    return series


def main():
    args = parse_args()
    window = args.velocity_window
    if window is None:
        window = 14 if args.dataset.startswith("sarscov2") else 60

    logfit = read_scaffolded(args.scaffolded)
    dates, by_date = read_frequencies(args.frequencies)

    variants = set(logfit) | {v for freqs in by_date.values() for v in freqs}
    aliasor = None
    if args.generation_time_pre_omicron is not None and args.variant_classification == "lineages":
        aliasor = load_aliasor(args.aliasing)
    tau_v = per_variant_tau(
        variants, args.generation_time, args.generation_time_pre_omicron,
        args.variant_classification, aliasor,
    )
    tau_bar = tau_bar_by_date(dates, by_date, tau_v, args.generation_time)

    location, variance = location_and_variance(dates, by_date, logfit)
    velocity = velocity_series(dates, location, tau_bar, window)
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
    mean_tau = sum(tau_bar[d] for d in dates) / len(dates) if dates else args.generation_time
    summary = {
        "dataset": args.dataset,
        "generation_time_days": args.generation_time,
        "mean_generation_time_days": mean_tau,
        "velocity_window_days": window,
        "avg_variance": float(np.mean(variance)),
        "avg_sd": float(np.mean([math.sqrt(v) for v in variance])),
        "avg_velocity": float(np.mean([v for _, v in velocity])) if velocity else None,
    }
    decimals = [ff_io.decimal_year(d) for d in dates]
    total_years = decimals[-1] - decimals[0]
    if total_years > 0:
        summary["simple_velocity"] = (
            (mean_tau / DAYS_PER_YEAR)
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
