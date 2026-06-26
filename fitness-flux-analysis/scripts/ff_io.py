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


def variant_growth_advantage_intervals(mlr, level=95):
    """variant -> (lower, upper) HDI bounds of the growth advantage ``ga``.

    Reads the ``HDI_{level}_lower`` / ``HDI_{level}_upper`` records of the ``ga``
    site for the primary location. The interval WIDTH (upper - lower) is a
    per-variant growth-rate uncertainty, available from posterior fits (NUTS,
    FullRank) but degenerate under a MAP point estimate. Variants without a ``ga``
    record (e.g. the pivot) are omitted.
    """
    location = primary_location(mlr)
    lo_key, hi_key = f"HDI_{level}_lower", f"HDI_{level}_upper"
    bounds = {}
    for record in mlr["data"]:
        if (
            record.get("value") is None
            or record.get("site") != "ga"
            or record.get("location") != location
        ):
            continue
        if record.get("ps") == lo_key:
            bounds.setdefault(record["variant"], {})["lo"] = record["value"]
        elif record.get("ps") == hi_key:
            bounds.setdefault(record["variant"], {})["hi"] = record["value"]
    return {v: (b["lo"], b["hi"]) for v, b in bounds.items() if "lo" in b and "hi" in b}


def variant_weekly_frequencies(mlr):
    """variant -> {date: weekly_raw_freq} for the primary-location series."""
    return _variant_date_series(mlr, site="weekly_raw_freq")


def variant_modeled_frequencies(mlr):
    """variant -> {date: MLR-modeled freq} (median) for the primary location."""
    return _variant_date_series(mlr, site="freq", ps="median")


def _variant_date_series(mlr, site, ps=None):
    location = primary_location(mlr)
    series = {}
    for record in mlr["data"]:
        if record.get("value") is None:
            continue
        if (
            record.get("site") == site
            and record.get("location") == location
            and (ps is None or record.get("ps") == ps)
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
