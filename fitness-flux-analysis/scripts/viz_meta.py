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


def parse_datasets(spec):
    """'id:Label,id2:Label2' -> [{"id","label"}, ...], order preserved."""
    out = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        dataset_id, _, label = item.partition(":")
        out.append({"id": dataset_id.strip(), "label": (label.strip() or dataset_id.strip())})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--datasets",
        default=None,
        help="override the dataset list as 'id:Label,...'; default is the full "
        "five-dataset manifest. Used by components that carry a subset (e.g. the "
        "similarity components, which carry only datasets with a clade-distance tree).",
    )
    parser.add_argument("--default", default=None, help="override the default dataset id")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    datasets = parse_datasets(args.datasets) if args.datasets else DATASETS
    default = args.default or (DEFAULT if not args.datasets else datasets[0]["id"])
    viz_io.write_json(args.output, {"datasets": datasets, "default": default}, indent=2)


if __name__ == "__main__":
    main()
