#!/usr/bin/env python3
"""Gather empirical variant frequencies and frequency-weighted mean dates.

Port of the "Gather frequencies" and "Mean date per variant" sections of
``fitness-flux.nb``. For each variant kept by the scaffolding step:

  * collect its weekly raw frequency across overlapping seasons, averaging the
    value where seasons overlap on a date, keeping dates above 0.005,
  * pad onto the union of all observed dates with zeros,
  * renormalize per date so the kept variants sum to 1 (``other`` excluded).

Outputs a long-format frequency table (positive entries only) and the
frequency-weighted mean decimal date per variant.
"""

import argparse
import csv

import ff_io


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="e.g. sarscov2_clades")
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument(
        "--scaffolded",
        required=True,
        help="scaffolded fitness TSV; its variants define the kept set",
    )
    parser.add_argument("--frequencies-output", required=True)
    parser.add_argument("--mean-date-output", required=True)
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.005,
        help="minimum averaged frequency to keep a (variant, date) entry",
    )
    return parser.parse_args()


def read_variants(scaffolded_path):
    with open(scaffolded_path) as handle:
        return [row["variant"] for row in csv.DictReader(handle, delimiter="\t")]


def gather(mlr_dir, dataset, variants, threshold):
    timepoints = ff_io.seasonal_timepoints(mlr_dir, dataset)

    # variant -> list of per-season {date: value} series
    per_season = {variant: [] for variant in variants}
    all_dates = set()
    for timepoint in timepoints:
        mlr = ff_io.load_mlr(mlr_dir, dataset, timepoint)
        season = ff_io.variant_weekly_frequencies(mlr)
        for variant in variants:
            series = season.get(variant, {})
            if len(series) > 1:  # notebook drops singleton series
                per_season[variant].append(series)
            all_dates.update(series)
    dates = sorted(all_dates)

    # average overlapping seasons per date, threshold, then pad to all dates
    freqs = {}
    for variant in variants:
        averaged = {}
        date_values = {}
        for series in per_season[variant]:
            for date, value in series.items():
                date_values.setdefault(date, []).append(value)
        for date, values in date_values.items():
            mean_value = sum(values) / len(values)
            if mean_value > threshold:
                averaged[date] = mean_value
        freqs[variant] = [averaged.get(date, 0.0) for date in dates]

    # renormalize per date so kept variants sum to 1
    for i in range(len(dates)):
        total = sum(freqs[variant][i] for variant in variants)
        if total > 0:
            for variant in variants:
                freqs[variant][i] /= total
    return dates, freqs


def mean_dates(dates, freqs, variants):
    decimal = [ff_io.decimal_year(date) for date in dates]
    result = {}
    for variant in variants:
        series = freqs[variant]
        total = sum(series)
        if total > 0:
            result[variant] = sum(d * f for d, f in zip(decimal, series)) / total
    return result


def main():
    args = parse_args()
    variants = read_variants(args.scaffolded)
    dates, freqs = gather(args.mlr_dir, args.dataset, variants, args.threshold)

    written = 0
    with open(args.frequencies_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["date", "variant", "frequency"])
        for i, date in enumerate(dates):
            for variant in variants:
                value = freqs[variant][i]
                if value > 0:
                    writer.writerow([date, variant, f"{value:.6f}"])
                    written += 1
    ff_io.log(f"Wrote {written} frequency rows to {args.frequencies_output}")

    means = mean_dates(dates, freqs, variants)
    with open(args.mean_date_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant", "mean_date"])
        for variant in variants:
            if variant in means:
                writer.writerow([variant, f"{round(means[variant], 4):g}"])
    ff_io.log(f"Wrote {len(means)} mean-date rows to {args.mean_date_output}")


if __name__ == "__main__":
    main()
