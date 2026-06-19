"""Shared helpers for the fitness-flux analysis scripts.

Ports the data-access conventions of the original ``fitness-flux.nb`` Mathematica
notebook: reading per-season MLR estimates, ordering seasonal timepoints, and
mapping season labels / calendar dates to decimal years.
"""

import json
import sys
from datetime import date
from glob import glob
from pathlib import Path


def log(message):
    print(message, file=sys.stderr)


def seasonal_timepoints(mlr_dir, dataset):
    """Return the sorted season suffixes available for ``dataset``.

    A dataset such as ``sarscov2_clades`` aggregates the per-season MLR runs
    ``mlr-estimates/sarscov2_clades_2020``, ``..._2020-21`` and so on. The
    notebook orders these by the directory listing, which is lexicographic;
    ``sorted`` reproduces that order.
    """
    prefix = f"{dataset}_"
    seasons = []
    for path in glob(str(Path(mlr_dir) / f"{prefix}*")):
        name = Path(path).name
        if (Path(path) / "mlr_results.json").exists():
            seasons.append(name[len(prefix):])
    return sorted(seasons)


def mlr_path(mlr_dir, dataset, timepoint):
    return Path(mlr_dir) / f"{dataset}_{timepoint}" / "mlr_results.json"


def load_mlr(mlr_dir, dataset, timepoint):
    with open(mlr_path(mlr_dir, dataset, timepoint)) as handle:
        return json.load(handle)


def primary_location(mlr):
    """The single non-hierarchical geographic location for a dataset.

    The notebook hardcoded ``"USA"``, but the upstream MLR runs now use other
    locations for some datasets (e.g. ``sarscov2_clades`` reports United
    Kingdom). Each dataset's metadata lists exactly one real location plus the
    pooled ``hierarchical`` entry, so we select the real one dynamically.
    """
    locations = [loc for loc in mlr["metadata"]["location"] if loc != "hierarchical"]
    if not locations:
        raise ValueError("No non-hierarchical location in MLR metadata")
    return locations[0]


def variant_growth_advantages(mlr):
    """variant -> growth advantage (linear ``ga``) for the median series.

    Mirrors ``variantFitnessesForTimePoint``: a variant present in the metadata
    but lacking a ``ga`` record defaults to 1. ``other`` is retained here and is
    excluded later where appropriate.
    """
    location = primary_location(mlr)
    variants = mlr["metadata"]["variants"]
    ga = {}
    for record in mlr["data"]:
        if record.get("value") is None:
            continue
        if (
            record.get("site") == "ga"
            and record.get("ps") == "median"
            and record.get("location") == location
        ):
            ga[record["variant"]] = record["value"]
    return {variant: ga.get(variant, 1.0) for variant in variants}


def variant_weekly_frequencies(mlr):
    """variant -> {date: weekly_raw_freq} for the primary-location series."""
    location = primary_location(mlr)
    series = {}
    for record in mlr["data"]:
        if record.get("value") is None:
            continue
        if (
            record.get("site") == "weekly_raw_freq"
            and record.get("location") == location
        ):
            series.setdefault(record["variant"], {})[record["date"]] = record["value"]
    return series


def timepoint_to_numeric(timepoint):
    """Season label -> midpoint decimal year (``"2020-21"`` -> ``2020.5``)."""
    parts = []
    for piece in timepoint.split("-"):
        value = int(piece)
        parts.append(value if value > 2000 else value + 2000)
    return sum(parts) / len(parts)


def decimal_year(date_str):
    """Calendar date string (``YYYY-MM-DD``) -> decimal year."""
    d = date.fromisoformat(date_str)
    year_start = date(d.year, 1, 1)
    next_year_start = date(d.year + 1, 1, 1)
    return d.year + (d - year_start).days / (next_year_start - year_start).days
