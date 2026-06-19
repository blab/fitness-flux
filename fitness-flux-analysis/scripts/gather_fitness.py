#!/usr/bin/env python3
"""Gather direct (per-season) variant fitness from MLR estimates.

Port of the "Gather fitnesses" section of ``fitness-flux.nb``. For each season
of a dataset (e.g. ``sarscov2_clades``) this reads the median USA growth
advantage (``ga``) per variant and records its natural log as the direct log
fitness, tagged with the season midpoint year. The ``other`` category is
dropped. Each season's values are relative to that season's own pivot; the
scaffolding step later stitches the seasons onto a common scale.
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
    parser.add_argument("--output", required=True, help="output TSV path")
    return parser.parse_args()


def gather_direct_fitness(mlr_dir, dataset):
    import math

    rows = []
    for timepoint in ff_io.seasonal_timepoints(mlr_dir, dataset):
        numeric = ff_io.timepoint_to_numeric(timepoint)
        mlr = ff_io.load_mlr(mlr_dir, dataset, timepoint)
        for variant, ga in ff_io.variant_growth_advantages(mlr).items():
            if variant == "other":
                continue
            rows.append((numeric, variant, math.log(ga)))
    return rows


def main():
    args = parse_args()
    rows = gather_direct_fitness(args.mlr_dir, args.dataset)
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["timepoint", "variant", "log_fitness"])
        for numeric, variant, log_fitness in rows:
            writer.writerow([f"{numeric:g}", variant, f"{log_fitness:.6f}"])
    ff_io.log(f"Wrote {len(rows)} direct-fitness rows to {args.output}")


if __name__ == "__main__":
    main()
