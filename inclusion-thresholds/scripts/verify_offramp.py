#!/usr/bin/env python3
"""Tune the collapse-time core off-ramp to maximize the lineage mutation->fitness
correlation (lineages are used for phylogenetic contrasts, not flux, so the
recent-window over-dissolution into 'other' is an accepted cost).

For each core floor (and a no-off-ramp baseline), rebuild the deltas pipeline at
collapse_threshold=500 across all 12 seasons with fast MAP fits, then report
Pearson/Spearman of the per-branch spike/s1/rbd mutation delta vs delta_log_fitness
(as in lineage-deltas-analysis/scripts/predictor_correlations.py). Pick the floor
that maximizes Pearson.
"""

import os
import subprocess
import sys

import pandas as pd
from scipy import stats

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "fitness-flux-analysis", "scripts"))
sys.path.insert(0, os.path.join(REPO, "lineage-deltas-analysis", "scripts"))
import ld_io  # noqa: E402

ALIAS = os.path.join(REPO, "inclusion-thresholds", "scratch", "alias_key.json")
MAP_CONFIG = os.path.join(REPO, "inclusion-thresholds", "scripts", "mlr-config-lineage-map.yaml")
MUT = os.path.join(REPO, "mutation-counts", "results", "sarscov2_lineages_mut_counts.tsv")
ROOT = os.path.join(REPO, "inclusion-thresholds", "scratch", "offramp_verify")
SEASONS = ["2020", "2020-21", "2021", "2021-22", "2022", "2022-23",
           "2023", "2023-24", "2024", "2024-25", "2025", "2025-26"]
THRESHOLD = "500"
FLOORS = [None, 50, 100, 200, 300, 500]  # None = baseline (no off-ramp)


def run(cmd):
    subprocess.run(cmd, check=True, cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def mode_name(floor):
    return "base" if floor is None else f"off_{floor}"


def build_fits(floor):
    mode = mode_name(floor)
    sdir = os.path.join(ROOT, mode, "seq")
    mdir = os.path.join(ROOT, mode, "mlr")
    for tp in SEASONS:
        ds = f"sarscov2_lineages_{tp}"
        prep = os.path.join(REPO, "sequence-counts", ds, "prepared_seq_counts.tsv")
        if not os.path.exists(prep):
            continue
        seq_out = os.path.join(sdir, ds)
        mlr_out = os.path.join(mdir, ds)
        os.makedirs(seq_out, exist_ok=True)
        os.makedirs(mlr_out, exist_ok=True)
        collapsed = os.path.join(seq_out, "collapsed_seq_counts.tsv")
        rel = os.path.join(seq_out, "variant_relationships.tsv")
        if os.path.exists(os.path.join(mlr_out, "mlr_results.json")):
            continue
        cmd = ["python", "scripts/collapse-lineage-counts.py", "--seq-counts", prep,
               "--collapse-threshold", THRESHOLD, "--aliasing", ALIAS,
               "--output-seq-counts", collapsed]
        if floor is not None:
            cmd += ["--min-core-count", str(floor)]
        run(cmd)
        run(["python", "scripts/prepare-pango-relationships.py", "--seq-counts", collapsed,
             "--output-relationships", rel])
        run(["python", "scripts/run-mlr-model.py", "--config", MAP_CONFIG, "--seq-path", collapsed,
             "--export-path", mlr_out, "--generation-time", "3.2", "--data-name", "mlr"])
    print(f"  {mode} fits ready", flush=True)
    return mdir, sdir


def branch_table(mdir, sdir):
    mut = pd.read_csv(MUT, sep="\t").set_index("lineage")
    rows = []
    for br in ld_io.branches(mdir, sdir):
        c, p = br["child"], br["parent"]
        if c in mut.index and p in mut.index:
            rows.append({
                "delta_log_fitness": br["delta_log_fitness"],
                "delta_spike": mut.loc[c, "spike_muts"] - mut.loc[p, "spike_muts"],
                "delta_s1": mut.loc[c, "s1_muts"] - mut.loc[p, "s1_muts"],
                "delta_rbd": mut.loc[c, "rbd_muts"] - mut.loc[p, "rbd_muts"],
            })
    return pd.DataFrame(rows)


def main():
    print(f"{'floor':>6} {'n':>5} {'spike_r':>8} {'s1_r':>7} {'rbd_r':>7} {'spike_rho':>10}", flush=True)
    for floor in FLOORS:
        mdir, sdir = build_fits(floor)
        df = branch_table(mdir, sdir)
        def corr(col):
            m = df[col].notna() & df["delta_log_fitness"].notna()
            return (stats.pearsonr(df.loc[m, col], df.loc[m, "delta_log_fitness"])[0],
                    stats.spearmanr(df.loc[m, col], df.loc[m, "delta_log_fitness"])[0], int(m.sum()))
        sp, s1, rb = corr("delta_spike"), corr("delta_s1"), corr("delta_rbd")
        label = "base" if floor is None else str(floor)
        print(f"{label:>6} {sp[2]:>5} {sp[0]:>8.3f} {s1[0]:>7.3f} {rb[0]:>7.3f} {sp[1]:>10.3f}", flush=True)


if __name__ == "__main__":
    main()
