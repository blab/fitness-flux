"""Shared helpers for the lineage-deltas analysis scripts.

Ports the data-access conventions of ``lineage-deltas.nb``: per-season variant
log fitness from MLR, per-season parent->child lineage relationships, and the
generic "branch" construction that pairs a change along a parent->child edge
with the corresponding change in (per-season) log fitness.
"""

import csv
import json
import sys
from glob import glob
from pathlib import Path

DATASET = "sarscov2_lineages"


def log(message):
    print(message, file=sys.stderr)


def seasonal_timepoints(mlr_dir, dataset=DATASET):
    prefix = f"{dataset}_"
    seasons = []
    for path in glob(str(Path(mlr_dir) / f"{prefix}*")):
        if (Path(path) / "mlr_results.json").exists():
            seasons.append(Path(path).name[len(prefix):])
    return sorted(seasons)


def _primary_location(metadata):
    locations = [loc for loc in metadata["location"] if loc != "hierarchical"]
    if not locations:
        raise ValueError("No non-hierarchical location in MLR metadata")
    return locations[0]


def variant_log_fitness(mlr_dir, timepoint, dataset=DATASET):
    """variant -> log(growth advantage) for one season (``other`` dropped)."""
    import math

    with open(Path(mlr_dir) / f"{dataset}_{timepoint}" / "mlr_results.json") as handle:
        mlr = json.load(handle)
    location = _primary_location(mlr["metadata"])
    result = {}
    for record in mlr["data"]:
        if record.get("value") is None:
            continue
        if (
            record.get("site") == "ga"
            and record.get("ps") == "median"
            and record.get("location") == location
            and record.get("variant") != "other"
        ):
            result[record["variant"]] = math.log(record["value"])
    return result


def variant_parents(seqcounts_dir, timepoint, dataset=DATASET):
    """child variant -> closest parent for one season (rootless edges dropped)."""
    path = (
        Path(seqcounts_dir)
        / f"{dataset}_{timepoint}"
        / "variant_relationships.tsv"
    )
    parents = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            parent = row.get("closest_parent") or ""
            if parent.strip():
                parents[row["variant"]] = parent
    return parents


def branches(mlr_dir, seqcounts_dir, dataset=DATASET):
    """Yield every parent->child edge that has log fitness on both endpoints.

    Each branch is a dict with the season, parent, child and the per-season
    change in log fitness. Callers attach their own count/predictor deltas by
    joining on parent and child.
    """
    rows = []
    for timepoint in seasonal_timepoints(mlr_dir, dataset):
        fitness = variant_log_fitness(mlr_dir, timepoint, dataset)
        parents = variant_parents(seqcounts_dir, timepoint, dataset)
        for child, parent in parents.items():
            if child in fitness and parent in fitness:
                rows.append(
                    {
                        "timepoint": timepoint,
                        "parent": parent,
                        "child": child,
                        "delta_log_fitness": fitness[child] - fitness[parent],
                    }
                )
    return rows


def read_tsv_map(path, key_col, value_col, delimiter="\t", cast=float):
    """Read a two-column lookup, skipping rows that fail to cast."""
    result = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter=delimiter):
            try:
                result[row[key_col]] = cast(row[value_col])
            except (ValueError, KeyError, TypeError):
                continue
    return result
