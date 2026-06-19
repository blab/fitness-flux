#!/usr/bin/env python3
"""
Build the frequency-panels component's data.json: per-season empirical (dots)
and modeled (lines) variant frequencies, bundled with the variant color table.
Empty cells become null; dates stay ISO strings.
"""
import argparse

import viz_io


def build_frequency_panels(seasonal_path):
    return [
        {
            "timepoint": row["timepoint"],
            "date": row["date"],
            "variant": row["variant"],
            "empirical": viz_io.num(row["empirical"]),
            "modeled": viz_io.num(row["modeled"]),
        }
        for row in viz_io.read_tsv(seasonal_path)
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seasonal", required=True)
    parser.add_argument("--colors", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    viz_io.write_json(
        args.output,
        {
            "seasonal": build_frequency_panels(args.seasonal),
            "colors": viz_io.build_colors(args.colors),
        },
    )


if __name__ == "__main__":
    main()
