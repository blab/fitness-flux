#!/usr/bin/env python3
"""Stitch per-season variant fitnesses into one scaffolded fitness scale.

Port of the "Update fitnesses" section of ``fitness-flux.nb``. Seasonal MLR
growth advantages (``ga``) are only defined relative to each season's pivot, so
overlapping variants are used to rescale successive seasons onto a common
scale. Working on the linear ``ga`` scale:

  scaffold = fold over seasons, starting from season 1's ga map; for each later
  season, divide its values by the mean ratio (newSeason[v] / accumulated[v])
  taken over shared variants (excluding ``other``), then merge into the
  accumulator by averaging the values of any variant present in both.

The final scaffolded fitness is ``log(ga)`` for each variant (``other``
dropped), rounded to 4 decimal places to match the notebook export.
"""

import argparse
import csv
import math
from statistics import mean

import ff_io


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="e.g. sarscov2_clades")
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument("--output", required=True, help="output TSV path")
    parser.add_argument(
        "--validate",
        help="optional reference TSV (variant\\tlog_fitness) to compare against",
    )
    return parser.parse_args()


def update_fitnesses(comparison, to_update):
    """Rescale ``to_update`` onto ``comparison`` via the mean shared-variant ratio."""
    shared = [v for v in (set(comparison) & set(to_update)) if v != "other"]
    ratios = [to_update[v] / comparison[v] for v in shared]
    mean_ratio = mean(ratios)
    return {v: value / mean_ratio for v, value in to_update.items()}


def combine_fitnesses(mapping1, mapping2):
    """Merge two variant->ga maps, averaging values for shared variants."""
    combined = {}
    for variant in set(mapping1) | set(mapping2):
        values = [m[variant] for m in (mapping1, mapping2) if variant in m]
        combined[variant] = mean(values)
    return combined


def scaffold(mlr_dir, dataset):
    timepoints = ff_io.seasonal_timepoints(mlr_dir, dataset)
    if not timepoints:
        raise SystemExit(f"No seasonal MLR estimates found for {dataset}")

    accumulated = ff_io.variant_growth_advantages(
        ff_io.load_mlr(mlr_dir, dataset, timepoints[0])
    )
    for timepoint in timepoints[1:]:
        new_season = ff_io.variant_growth_advantages(
            ff_io.load_mlr(mlr_dir, dataset, timepoint)
        )
        rescaled = update_fitnesses(accumulated, new_season)
        accumulated = combine_fitnesses(accumulated, rescaled)

    return {
        variant: math.log(ga)
        for variant, ga in accumulated.items()
        if variant != "other"
    }


def read_reference(path):
    reference = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reference[row["variant"]] = float(row["log_fitness"])
    return reference


def validate(scaffolded, reference_path):
    reference = read_reference(reference_path)
    keys = set(scaffolded) | set(reference)
    worst = 0.0
    worst_variant = None
    missing = []
    for variant in sorted(keys):
        if variant not in scaffolded or variant not in reference:
            missing.append(variant)
            continue
        diff = abs(round(scaffolded[variant], 4) - reference[variant])
        if diff > worst:
            worst, worst_variant = diff, variant
    ff_io.log(
        f"Validation vs {reference_path}: max |Δlog_fitness| = {worst:.4f} "
        f"(variant {worst_variant}); {len(missing)} unmatched variants"
    )
    if missing:
        ff_io.log(f"  unmatched: {', '.join(missing)}")
    return worst, missing


def main():
    args = parse_args()
    scaffolded = scaffold(args.mlr_dir, args.dataset)
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant", "log_fitness"])
        for variant, log_fitness in scaffolded.items():
            writer.writerow([variant, f"{round(log_fitness, 4):g}"])
    ff_io.log(f"Wrote {len(scaffolded)} scaffolded-fitness rows to {args.output}")
    if args.validate:
        validate(scaffolded, args.validate)


if __name__ == "__main__":
    main()
