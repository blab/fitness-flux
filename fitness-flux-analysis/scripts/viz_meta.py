#!/usr/bin/env python3
"""
Emit a component's meta.json manifest: the datasets it carries (id + human
label), in display order, plus the default dataset. Consumed by the dashboard
selector and the dev harness so hosts don't hardcode the dataset list (static
hosts can't list a directory).

The dataset ids follow the analysis pipeline naming; the labels and order are a
presentation concern curated here. Both components currently carry the same
datasets, so this manifest is shared.
"""
import argparse

import viz_io

DATASETS = [
    {"id": "sarscov2_clades", "label": "SARS-CoV-2 clades"},
    {"id": "sarscov2_lineages", "label": "SARS-CoV-2 lineages"},
    {"id": "h3n2_clades", "label": "H3N2 clades"},
]
DEFAULT = "sarscov2_clades"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    viz_io.write_json(args.output, {"datasets": DATASETS, "default": DEFAULT}, indent=2)


if __name__ == "__main__":
    main()
