#!/usr/bin/env python3
"""Join per-variant mutation counts with scaffolded fitness and mean date.

Port of the "Plotting mutations and fitness" inputs of ``fitness-flux.nb``.
Mutation counts (receptor-binding protein amino-acid substitutions per clade,
from ``source-data/{dataset}_mut_counts.tsv``) are joined to the scaffolded log
fitness and frequency-weighted mean date so the dashboard can show mutation
count versus fitness and mutation accumulation through time.
"""

import argparse
import csv

import ff_io


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--mut-counts", required=True)
    parser.add_argument("--scaffolded", required=True)
    parser.add_argument("--mean-date", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_map(path, key, value, cast=float):
    result = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                result[row[key]] = cast(row[value])
            except (ValueError, KeyError):
                continue
    return result


def main():
    args = parse_args()
    mut_counts = read_map(args.mut_counts, "clade", "mutations")
    log_fitness = read_map(args.scaffolded, "variant", "log_fitness")
    mean_date = read_map(args.mean_date, "variant", "mean_date")

    rows = []
    for variant, fitness in log_fitness.items():
        rows.append(
            {
                "variant": variant,
                "mut_count": mut_counts.get(variant),
                "log_fitness": fitness,
                "mean_date": mean_date.get(variant),
            }
        )

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["variant", "mut_count", "log_fitness", "mean_date"])
        for row in rows:
            writer.writerow(
                [
                    row["variant"],
                    "" if row["mut_count"] is None else f"{row['mut_count']:g}",
                    f"{row['log_fitness']:g}",
                    "" if row["mean_date"] is None else f"{row['mean_date']:g}",
                ]
            )
    ff_io.log(f"Wrote {len(rows)} mutation-fitness rows to {args.output}")


if __name__ == "__main__":
    main()
