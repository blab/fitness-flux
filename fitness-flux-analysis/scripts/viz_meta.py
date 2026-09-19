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
    {"id": "h1n1pdm_clades", "label": "H1N1pdm clades"},
    {"id": "vic_clades", "label": "Vic clades"},
]
DEFAULT = "sarscov2_clades"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--include",
        default=None,
        help="comma-separated dataset ids to keep (in the given order); default keeps all",
    )
    args = parser.parse_args()

    datasets = DATASETS
    default = DEFAULT
    if args.include:
        wanted = [d.strip() for d in args.include.split(",") if d.strip()]
        by_id = {d["id"]: d for d in DATASETS}
        datasets = [by_id[i] for i in wanted if i in by_id]
        if datasets:
            default = datasets[0] if DEFAULT not in {d["id"] for d in datasets} else DEFAULT

    viz_io.write_json(args.output, {"datasets": datasets, "default": default}, indent=2)


if __name__ == "__main__":
    main()
