#!/usr/bin/env python3
"""
Build the fitness-flux component's data.json: each empirical-frequency point
above a threshold whose variant has a scaffolded log fitness, joined to that
fitness and bundled with the variant color table. Emits per record
{variant, date, emp_freq, log_fitness} (emp_freq = empirical frequency,
log_fitness = scaffolded log fitness).

The join is frozen here so the browser component stays pure. Mirrors the inline
join in the original viz/fitness-flux.html renderFitnessFlux. Dates stay ISO
strings; the component constructs Date objects.
"""
import argparse

import viz_io


def build_fitness_flux(frequencies_path, scaffolded_path, min_freq):
    log_fit = {}
    for row in viz_io.read_tsv(scaffolded_path):
        fit = viz_io.num(row["log_fitness"])
        if fit is not None:
            log_fit[row["variant"]] = fit
    flux = []
    for row in viz_io.read_tsv(frequencies_path):
        variant = row["variant"]
        freq = viz_io.num(row["frequency"])
        if freq is None or freq <= min_freq or variant not in log_fit:
            continue
        flux.append(
            {
                "variant": variant,
                "date": row["date"],
                "emp_freq": freq,
                "log_fitness": log_fit[variant],
            }
        )
    return flux


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frequencies", required=True)
    parser.add_argument("--scaffolded", required=True)
    parser.add_argument("--colors", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-freq", type=float, default=0.01)
    args = parser.parse_args()

    viz_io.write_json(
        args.output,
        {
            "fitnessFlux": build_fitness_flux(args.frequencies, args.scaffolded, args.min_freq),
            "colors": viz_io.build_colors(args.colors),
        },
    )


if __name__ == "__main__":
    main()
