#!/usr/bin/env python3
"""Stitch per-window variant fitnesses into one scaffolded fitness scale.

Each sliding-window MLR fit estimates a variant's fitness ``f = log(ga)`` only
relative to that window's pivot, so every window carries its own arbitrary
additive zero. Scaffolding recovers a single fitness per variant by treating
this as a weighted two-way additive model: each per-window estimate is

    f_{i,w} ~= f_i - c_w,

a variant effect ``f_i`` (its global fitness) minus a window effect ``c_w``
(window w's offset). We choose the ``f_i`` and ``c_w`` that jointly minimize the
abundance-weighted squared error over every window,

    minimize  sum_{i,w} a_{i,w} (f_{i,w} - f_i + c_w)^2,

weighting each estimate by ``a_{i,w}``, the area under variant i's modeled-
frequency curve in window w (AUC) — so a window where a variant is a near-
extinct relic, its MLR estimate poorly constrained, contributes negligibly.

The optimum is two interleaved abundance-weighted means: each variant's fitness
is the weighted mean of its offset-corrected estimates over the windows it
appears in, and each window's offset is the weighted-mean gap between the global
scale and that window's estimates. We solve by alternating the two to
convergence (the overlap between windows ties them into one connected scale, so
the relative fitnesses are determined up to a single global constant). Finally
the scale is shifted so the founding clade — the least-fit, our baseline — sits
at 0, leaving each variant's value as its cumulative fitness flux relative to it
(``f_i - f_0 >= 0``). Rounded to 4 decimals.
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


def window_fitnesses(mlr):
    """variant -> fitness ``f = log(ga)`` for the window (relative to its pivot),
    excluding the pooled ``other`` category, whose aggregate fitness is not a
    meaningful per-variant quantity and must not anchor a window's offset."""
    return {
        variant: math.log(ga)
        for variant, ga in ff_io.variant_growth_advantages(mlr).items()
        if variant != "other"
    }


def _weighted_mean(pairs):
    """Mean of ``(value, weight)`` pairs; falls back to the unweighted mean when
    no entry carries weight (e.g. a variant present only in relic windows)."""
    total = sum(w for _, w in pairs)
    if total > 0:
        return sum(value * w for value, w in pairs) / total
    return mean([value for value, _ in pairs])


def two_way_fit(windows, max_iter=10000, tol=1e-12):
    """Jointly fit a global fitness ``f_i`` per variant and an offset ``c_w`` per
    window to all per-window estimates by abundance-weighted least squares,
    minimizing ``sum_{i,w} a_{i,w} (f_{i,w} - f_i + c_w)^2``.

    ``windows`` is a list of ``(fitness, weight)`` dicts, one per window, mapping
    variant -> ``f_{i,w}`` and variant -> ``a_{i,w}`` (AUC). The optimum is the
    fixed point of two interleaved weighted means, reached by alternating them:
    each window offset is the weighted-mean gap to the current global scale, and
    each variant fitness is the weighted mean of its offset-corrected estimates.
    The absolute level is gauge-free (anchored later); iteration stops once no
    fitness moves by more than ``tol``.
    """
    # variant -> list of (window index, f_{i,w}, a_{i,w}); and the per-window dual
    by_variant = {}
    for w, (fitness_w, weight_w) in enumerate(windows):
        for variant, f in fitness_w.items():
            by_variant.setdefault(variant, []).append((w, f, weight_w.get(variant, 0.0)))
    by_window = [
        [(variant, f, weight_w.get(variant, 0.0)) for variant, f in fitness_w.items()]
        for fitness_w, weight_w in windows
    ]

    offsets = [0.0] * len(windows)
    fitness = {
        variant: _weighted_mean([(f, a) for _, f, a in entries])
        for variant, entries in by_variant.items()
    }
    for _ in range(max_iter):
        for w, entries in enumerate(by_window):
            offsets[w] = _weighted_mean([(fitness[v] - f, a) for v, f, a in entries])
        max_delta = 0.0
        for variant, entries in by_variant.items():
            updated = _weighted_mean([(f + offsets[w], a) for w, f, a in entries])
            max_delta = max(max_delta, abs(updated - fitness[variant]))
            fitness[variant] = updated
        if max_delta < tol:
            break
    return fitness


def scaffold(mlr_dir, dataset):
    timepoints = ff_io.seasonal_timepoints(mlr_dir, dataset)
    if not timepoints:
        raise SystemExit(f"No seasonal MLR estimates found for {dataset}")

    windows = []
    for timepoint in timepoints:
        mlr = ff_io.load_mlr(mlr_dir, dataset, timepoint)
        windows.append((window_fitnesses(mlr), window_weights(mlr)))

    fitness = two_way_fit(windows)
    # Anchor the founding clade (the least-fit baseline) to 0, so each value is the
    # cumulative fitness flux relative to it (f_i - f_0 >= 0).
    floor = min(fitness.values())
    return {variant: value - floor for variant, value in fitness.items()}


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
