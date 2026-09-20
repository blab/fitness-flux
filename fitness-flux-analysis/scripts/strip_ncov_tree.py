#!/usr/bin/env python3
"""Strip the heavy ``tree`` from a Nextstrain auspice JSON, keeping ``meta`` (and ``version``).

``fitness_flux_colors`` only reads ``meta.colorings`` (the ``clade_membership`` scale), so the
tree-stripped metadata is committed under ``source-data/`` as a small pinned input instead of the
full multi-megabyte auspice tree. See ``rules/fitness_flux_analysis.smk`` for how it is refreshed.
"""
import argparse
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, help="full Nextstrain auspice tree JSON")
    parser.add_argument("--output", required=True, help="output path for tree-stripped metadata JSON")
    args = parser.parse_args()

    with open(args.tree) as handle:
        auspice = json.load(handle)
    auspice.pop("tree", None)
    with open(args.output, "w") as handle:
        json.dump(auspice, handle)


if __name__ == "__main__":
    main()
