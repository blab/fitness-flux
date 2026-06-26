#!/usr/bin/env python3
"""Sweep `clade_min_seq` for representative MLR datasets and fit each with MAP.

For each (dataset, threshold theta), reuse the cached, threshold-independent
`sequence-counts/{dataset}/seq_counts.tsv`, re-run `prepare-data.py` (date-window
+ collapse clades below theta into "other"; the pivot is force-included so it
survives), and fit a fast MAP MLR into
`inclusion-thresholds/scratch/{dataset}/thresh_{theta}/`. Lineage datasets
instead vary `collapse_threshold` (hierarchical roll-up into the parent lineage).

This does NOT touch the committed config or the canonical sequence-counts/ and
mlr-estimates/ outputs. Run from the repo root: idempotent (skips a (dataset,
theta) whose mlr_results.json already exists), so it can be resumed/batched.
"""

import argparse
import os
import subprocess

import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
ANALYSIS = os.path.join(REPO, "inclusion-thresholds")
SCRATCH = os.path.join(ANALYSIS, "scratch")
MAP_CONFIG = os.path.join(ANALYSIS, "scripts", "mlr-config-map.yaml")
LINEAGE_MAP_CONFIG = os.path.join(ANALYSIS, "scripts", "mlr-config-lineage-map.yaml")
ALIAS_FILE = os.path.join(SCRATCH, "alias_key.json")

# Representative windows spanning sequence volume. `current` is the production
# threshold (clade_min_seq for clades, collapse_threshold for lineages), used for
# reference markers. Clade datasets collapse rare clades into "other" (MAP fits);
# lineage datasets roll rare Pango lineages up into their PARENT via
# collapse-lineage-counts.py and use fast MAP fits (pivot auto-selected, =None) for
# growth advantages + frequencies; the noise metric is computed analytically from
# the counts (see score_lineage_sweep.py), so a point estimate suffices.
DATASETS = {
    "sarscov2_clades_2020":      dict(kind="clade",   min_date="2020-01-01", max_date="2020-12-31", pivot="19B", gen=3.2, rfw=7,  current=1000),
    "sarscov2_clades_2021-22":   dict(kind="clade",   min_date="2021-07-01", max_date="2022-06-30", pivot="21K", gen=3.2, rfw=7,  current=1000),
    "sarscov2_clades_2024":      dict(kind="clade",   min_date="2024-01-01", max_date="2024-12-31", pivot="24A", gen=3.2, rfw=7,  current=500),
    "h3n2_clades_2023-24":       dict(kind="clade",   min_date="2023-01-01", max_date="2024-12-31", pivot="J.2", gen=3.2, rfw=14, current=100),
    "sarscov2_lineages_2020-21": dict(kind="lineage", min_date="2020-07-01", max_date="2021-06-30", pivot=None, gen=3.2, rfw=7,  current=1000),
    "sarscov2_lineages_2021-22": dict(kind="lineage", min_date="2021-07-01", max_date="2022-06-30", pivot=None, gen=3.2, rfw=7,  current=1000),
    "sarscov2_lineages_2023":    dict(kind="lineage", min_date="2023-01-01", max_date="2023-12-31", pivot=None, gen=3.2, rfw=7,  current=1000),
    "sarscov2_lineages_2025":    dict(kind="lineage", min_date="2025-01-01", max_date="2025-12-31", pivot=None, gen=3.2, rfw=7,  current=1000),
}

COUNT_GRID = [50, 100, 200, 500, 1000, 2000, 5000]
FREQ_PCTS = dict(clade=[0.05, 0.1, 0.2, 0.5, 1.0, 2.0],
                 lineage=[0.02, 0.05, 0.1, 0.2, 0.5, 1.0])


def windowed_N(name):
    """Total sequences in the window (collapse-invariant) from the canonical output."""
    f = os.path.join(REPO, "sequence-counts", name, "collapsed_seq_counts.tsv")
    return int(pd.read_csv(f, sep="\t")["sequences"].sum())


def thetas_for(kind, N, current):
    grid = {t for t in COUNT_GRID if t <= 0.5 * N}
    grid |= {max(1, round(p / 100 * N)) for p in FREQ_PCTS[kind]}
    grid.add(current)
    return sorted(grid)


def run(cmd, log):
    print("  +", " ".join(str(c) for c in cmd), flush=True)
    with open(log, "w") as lh:
        subprocess.run(cmd, check=True, stdout=lh, stderr=subprocess.STDOUT, cwd=REPO)


def sweep_one(ds, theta):
    name = ds["name"]
    outdir = os.path.join(SCRATCH, name, f"thresh_{theta}")
    if os.path.exists(os.path.join(outdir, "mlr_results.json")):
        print(f"  skip {name} theta={theta} (exists)", flush=True)
        return
    os.makedirs(outdir, exist_ok=True)
    seq_counts = os.path.join(REPO, "sequence-counts", name, "seq_counts.tsv")
    prepared = os.path.join(outdir, "prepared_seq_counts.tsv")

    cmd = ["python", "scripts/prepare-data.py",
           "--seq-counts", seq_counts,
           "--min-date", ds["min_date"], "--max-date", ds["max_date"],
           "--output-seq-counts", prepared]
    if ds["kind"] == "clade":
        cmd += ["--clade-min-count", str(theta), "--force-include-clades", ds["pivot"]]
    else:
        cmd += ["--clade-min-count", "1"]
    run(cmd, os.path.join(outdir, "prepare.log"))

    mlr_input = prepared
    if ds["kind"] == "lineage":
        collapsed = os.path.join(outdir, "collapsed_seq_counts.tsv")
        run(["python", "scripts/collapse-lineage-counts.py",
             "--seq-counts", prepared, "--collapse-threshold", str(theta),
             "--aliasing", ALIAS_FILE, "--output-seq-counts", collapsed],
            os.path.join(outdir, "collapse.log"))
        mlr_input = collapsed

    config = LINEAGE_MAP_CONFIG if ds["kind"] == "lineage" else MAP_CONFIG
    cmd = ["python", "scripts/run-mlr-model.py", "--config", config,
           "--seq-path", mlr_input, "--export-path", outdir,
           "--generation-time", str(ds["gen"]),
           "--raw-freq-window", str(ds["rfw"]), "--data-name", "mlr"]
    if ds["pivot"]:  # clades pin a pivot; lineages auto-select (pivot=None)
        cmd += ["--pivot", ds["pivot"]]
    run(cmd, os.path.join(outdir, "mlr.log"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default="all", help="dataset name or 'all'")
    ap.add_argument("--theta", type=int, help="run a single threshold (smoke test)")
    args = ap.parse_args()

    names = list(DATASETS) if args.dataset == "all" else [args.dataset]
    for name in names:
        ds = {**DATASETS[name], "name": name}
        N = windowed_N(name)
        thetas = [args.theta] if args.theta else thetas_for(ds["kind"], N, ds["current"])
        print(f"=== {name}: N={N:,}, thetas={thetas} ===", flush=True)
        for theta in thetas:
            sweep_one(ds, theta)


if __name__ == "__main__":
    main()
