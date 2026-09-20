#!/usr/bin/env python3
"""Export per-season empirical and MLR-modeled variant frequencies.

Drives the per-season small-multiple frequency panels (the notebook's
{dataset}_time_vs_frequency figure), which overlay empirical weekly frequencies
against the MLR-modeled frequencies for each variant in each season. This is the
raw per-season data, distinct from the cross-season gathered frequencies in
gather_frequencies.py.

Output {dataset}_seasonal_frequencies.tsv: timepoint, date, variant, empirical,
modeled (one row per variant-date where either series is present). The aggregated
"other" category is excluded and the named clades are renormalized per date so
they sum to 1 (empirical and modeled each renormalized by their own named total),
matching the named-clade-only convention used elsewhere; individual clades are
shown only above --min-peak.
"""

import argparse
import csv

import ff_io


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument(
        "--min-peak",
        type=float,
        default=0.05,
        help="keep a variant in a season only if its peak frequency reaches this",
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def renormalized(series_by_variant, named):
    """Per-date renormalization of a {variant: {date: value}} map over named clades.

    Returns {variant: {date: value/named_total}} using each date's sum over named
    clades (non-null values) as the denominator; dates whose named total is <= 0 are
    dropped. "other" is excluded from both the numerator and the denominator.
    """
    totals = {}
    for variant in named:
        for date, value in series_by_variant.get(variant, {}).items():
            if value is not None:
                totals[date] = totals.get(date, 0.0) + value
    out = {}
    for variant in named:
        renorm = {}
        for date, value in series_by_variant.get(variant, {}).items():
            total = totals.get(date, 0.0)
            if value is not None and total > 0:
                renorm[date] = value / total
        out[variant] = renorm
    return out


def main():
    args = parse_args()
    rows = []
    for timepoint in ff_io.seasonal_timepoints(args.mlr_dir, args.dataset):
        mlr = ff_io.load_mlr(args.mlr_dir, args.dataset, timepoint)
        # Named clades only (exclude the aggregated "other" bucket), each series
        # renormalized per date so the named clades sum to 1.
        named = [v for v in mlr["metadata"]["variants"] if v != "other"]
        empirical = renormalized(ff_io.variant_weekly_frequencies(mlr), named)
        modeled = renormalized(ff_io.variant_modeled_frequencies(mlr), named)
        for variant in named:
            emp, mod = empirical.get(variant, {}), modeled.get(variant, {})
            peak = max([v for v in list(emp.values()) + list(mod.values())], default=0.0)
            if peak < args.min_peak:
                continue
            for date in sorted(set(emp) | set(mod)):
                rows.append((timepoint, date, variant, emp.get(date), mod.get(date)))

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["timepoint", "date", "variant", "empirical", "modeled"])
        for timepoint, date, variant, emp, mod in rows:
            writer.writerow([
                timepoint,
                date,
                variant,
                "" if emp is None else f"{emp:.4f}",
                "" if mod is None else f"{mod:.4f}",
            ])
    ff_io.log(f"Wrote {len(rows)} seasonal-frequency rows to {args.output}")


if __name__ == "__main__":
    main()
