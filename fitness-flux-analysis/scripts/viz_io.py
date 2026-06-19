"""Shared helpers for building the portable viz components' JSON data.

Each figure's data.json is a complete, render-ready payload: the figure data
plus the variant color table (data.colors), so the browser component carries no
parsing, join, or color logic and stays pure.
"""
import csv
import json
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(value):
    return None if value is None or value == "" else float(value)


def build_colors(colors_path):
    return [
        {
            "variant": row["variant"],
            "color": row["color"],
            "display_name": row["display_name"],
            "is_major": str(row["is_major"]).strip().lower() == "true",
            "order": int(float(row["order"])),
        }
        for row in read_tsv(colors_path)
    ]


def write_json(path, obj, indent=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        if indent is None:
            json.dump(obj, handle, separators=(",", ":"))
        else:
            json.dump(obj, handle, indent=indent)
        handle.write("\n")
