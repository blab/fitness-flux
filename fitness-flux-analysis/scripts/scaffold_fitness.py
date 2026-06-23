#!/usr/bin/env python3
"""Stitch per-season variant fitnesses into one scaffolded fitness scale.

Port of the "Update fitnesses" section of ``fitness-flux.nb``. Seasonal MLR
growth advantages (``ga``) are only defined relative to each season's pivot, so
overlapping variants are used to rescale successive seasons onto a common
scale. Working on the linear ``ga`` scale:

  scaffold = fold over seasons, starting from season 1's ga map; for each later
  season, divide its values by the (abundance-weighted) mean ratio
  (newSeason[v] / accumulated[v]) over shared variants (excluding ``other``),
  then merge into the accumulator by an abundance-weighted average of any
  variant present in both.

Each season's contribution to a variant is weighted by how much that variant
circulated in the season — the area under its modeled-frequency curve (AUC) — so
seasons where a variant is a near-extinct relic (an unconstrained MLR estimate)
carry negligible weight. The scaffolded fitness is ``log(ga)`` for each variant
(``other`` dropped), then shifted so the founding clade — the least-fit baseline,
i.e. the minimum — sits at 0; each value is thus the cumulative fitness flux
relative to the founding variant (``f_i - f_0 >= 0``). Rounded to 4 decimals.
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


def window_weights(mlr):
    """variant -> area under its modeled-frequency curve over the window.

    A smooth measure (trapezoidal integral over decimal years, in freq*year
    units) of how much a variant circulated in the season, used to down-weight
    seasons where a variant is a near-extinct relic. Variants with fewer than two
    dated points integrate to 0.
    """
    weights = {}
    for variant, series in ff_io.variant_modeled_frequencies(mlr).items():
        points = sorted(
            (ff_io.decimal_year(date), freq) for date, freq in series.items()
        )
        weights[variant] = sum(
            0.5 * (points[i][1] + points[i + 1][1]) * (points[i + 1][0] - points[i][0])
            for i in range(len(points) - 1)
        )
    return weights


def update_fitnesses(comparison, to_update, weights):
    """Rescale ``to_update`` onto ``comparison`` via the abundance-weighted mean
    shared-variant ratio (abundant shared variants anchor the rescale; falls back
    to the plain mean if no shared variant carries weight)."""
    shared = [v for v in (set(comparison) & set(to_update)) if v != "other"]
    ratios = [to_update[v] / comparison[v] for v in shared]
    ws = [weights.get(v, 0.0) for v in shared]
    total = sum(ws)
    mean_ratio = sum(r * w for r, w in zip(ratios, ws)) / total if total > 0 else mean(ratios)
    return {v: value / mean_ratio for v, value in to_update.items()}


def combine_fitnesses(acc_ga, acc_w, new_ga, new_w):
    """Merge a new (rescaled) season into the accumulator, averaging each shared
    variant's ga weighted by its per-season frequency-area and tracking the
    accumulated weight."""
    combined_ga, combined_w = {}, {}
    for variant in set(acc_ga) | set(new_ga):
        ag, aw = acc_ga.get(variant), acc_w.get(variant, 0.0)
        ng, nw = new_ga.get(variant), new_w.get(variant, 0.0)
        if variant in acc_ga and variant in new_ga:
            total = aw + nw
            combined_ga[variant] = (ag * aw + ng * nw) / total if total > 0 else (ag + ng) / 2
            combined_w[variant] = total
        elif variant in acc_ga:
            combined_ga[variant], combined_w[variant] = ag, aw
        else:
            combined_ga[variant], combined_w[variant] = ng, nw
    return combined_ga, combined_w


def scaffold(mlr_dir, dataset):
    timepoints = ff_io.seasonal_timepoints(mlr_dir, dataset)
    if not timepoints:
        raise SystemExit(f"No seasonal MLR estimates found for {dataset}")

    first = ff_io.load_mlr(mlr_dir, dataset, timepoints[0])
    accumulated = ff_io.variant_growth_advantages(first)
    seed_w = window_weights(first)
    weight = {variant: seed_w.get(variant, 0.0) for variant in accumulated}
    for timepoint in timepoints[1:]:
        mlr = ff_io.load_mlr(mlr_dir, dataset, timepoint)
        new_ga = ff_io.variant_growth_advantages(mlr)
        new_w = window_weights(mlr)
        rescaled = update_fitnesses(accumulated, new_ga, new_w)
        accumulated, weight = combine_fitnesses(accumulated, weight, rescaled, new_w)

    scaffolded = {
        variant: math.log(ga)
        for variant, ga in accumulated.items()
        if variant != "other"
    }
    # Anchor the founding clade (the least-fit baseline) to 0, so each value is the
    # cumulative fitness flux relative to it (f_i - f_0 >= 0).
    floor = min(scaffolded.values())
    return {variant: value - floor for variant, value in scaffolded.items()}


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
