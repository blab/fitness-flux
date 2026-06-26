#!/usr/bin/env python3
"""Score the lineage `collapse_threshold` sweep.

Unlike the clade case (collapse into a universal "other"), a rare Pango lineage is
rolled UP into its parent — a coherent ancestor — so the two failure modes live on
different scales:

* NOISE (threshold too low): a retained lineage's growth-rate estimate is imprecise
  when it has few sequences. Quantified by a closed-form Fisher-information standard
  error from the counts, SE_v = 1/sqrt(sum_t N_t f_vt(1-f_vt)(t-tbar)^2) with t in
  generations -- the logistic-slope precision, which is governed by absolute COUNT
  (~1/sqrt(n)). This is deterministic (no fitted posterior, no ADVI convergence
  artifacts); NUTS at the recommended threshold confirms it tracks real HDI.
* OVER-COLLAPSE (threshold too high): merging a sub-lineage into its parent erases
  any real fitness difference between them. Measured from a fine-resolution
  reference fit by re-grouping fine lineages under the ancestor they would collapse
  to at theta and taking the frequency-weighted variance of their growth advantages
  ("fitness signal lost"). This is governed by FREQUENCY x fitness-difference, and
  is gentle because parent and child fitnesses are usually close.

Also reports the confounded COARSE cross-entropy (scored on each fit's own labels,
which rewards collapsing) and a FIXED-partition CE (every fit scored on the common
finest partition, a collapsed lineage predicted via its retained ancestor's
trajectory split uniformly among co-collapsed fine lineages) for continuity with
the clade analysis.

Writes results/lineage_sweep_metrics.tsv (per threshold) and
results/lineage_peaks.tsv (per retained lineage: count, freq, ga, se_proxy).
"""

import csv
import json
import math
import os
import sys
from datetime import date as _date

import pandas as pd
from pango_aliasor.aliasor import Aliasor

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRATCH = os.path.join(REPO, "inclusion-thresholds", "scratch")
RESDIR = os.path.join(REPO, "inclusion-thresholds", "results")
ALIAS_FILE = os.path.join(SCRATCH, "alias_key.json")
sys.path.insert(0, os.path.join(REPO, "fitness-flux-analysis", "scripts"))
sys.path.insert(0, os.path.dirname(__file__))

from ff_io import (  # noqa: E402
    variant_growth_advantages,
    variant_modeled_frequencies,
)
from sweep_clade_min_seq import DATASETS  # noqa: E402

EPS = 1e-6
ALIASOR = Aliasor(ALIAS_FILE)
_DEPTH = {}


def depth(lineage):
    if lineage == "other":
        return 0
    if lineage not in _DEPTH:
        _DEPTH[lineage] = len(ALIASOR.uncompress(lineage).split("."))
    return _DEPTH[lineage]


def collapse_grouping(fine_totals, theta):
    """Replicate collapse-lineage-counts.py: map each fine lineage to the ancestor
    it rolls up into at `theta` (rootless basal lineages -> "other")."""
    label = {l: l for l in fine_totals}          # fine lineage -> current label
    total = dict(fine_totals)                     # current label -> window total

    def low():
        return {v for v, c in total.items() if 0 < c < theta and v != "other"}

    cur = low()
    max_depth = max((depth(v) for v in cur), default=0)
    for d in range(max_depth, 0, -1):
        for v in [x for x in low() if depth(x) == d]:
            parent = ALIASOR.parent(v) or "other"
            total[parent] = total.get(parent, 0) + total.get(v, 0)
            total[v] = 0
            for fine, lab in label.items():
                if lab == v:
                    label[fine] = parent
    return label


def read_counts(run_dir):
    """(by_date {date:{variant:count}}, totals {variant:count}, N) from the MLR input."""
    path = os.path.join(run_dir, "collapsed_seq_counts.tsv")
    df = pd.read_csv(path, sep="\t")
    df["date"] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%d")
    df = df.groupby(["date", "variant"], as_index=False)["sequences"].sum()
    by_date = {}
    for date, sub in df.groupby("date"):
        by_date[date] = dict(zip(sub["variant"], sub["sequences"]))
    totals = df.groupby("variant")["sequences"].sum().astype(int).to_dict()
    return by_date, totals, int(df["sequences"].sum())


def se_proxies(by_date, min_date, gen):
    """Closed-form Fisher-information SE of each variant's logistic growth rate,
    from the daily counts: SE_v = 1/sqrt(sum_t N_t f_vt(1-f_vt)(t-tbar)^2), t in
    generations. Precision is set by absolute count (and temporal spread)."""
    d0 = _date.fromisoformat(min_date)
    acc = {}  # variant -> [sum_w, sum_w_t, sum_w_t2]
    for ds, counts in by_date.items():
        Nd = sum(counts.values())
        if Nd <= 0:
            continue
        t = (_date.fromisoformat(ds) - d0).days / gen
        for v, n in counts.items():
            f = n / Nd
            w = Nd * f * (1 - f)
            if w <= 0:
                continue
            a = acc.setdefault(v, [0.0, 0.0, 0.0])
            a[0] += w
            a[1] += w * t
            a[2] += w * t * t
    out = {}
    for v, (sw, swt, swt2) in acc.items():
        info = swt2 - swt * swt / sw if sw > 0 else 0.0
        out[v] = 1.0 / math.sqrt(info) if info > 1e-12 else float("inf")
    return out


def coarse_ce(by_date, totals, N, modeled):
    ll = 0.0
    for date, counts in by_date.items():
        for v, n in counts.items():
            f = max(modeled.get(v, {}).get(date, EPS), EPS)
            ll += n * math.log(f)
    return -ll / N


def fixed_ce(fine_by_date, N, modeled, grouping):
    """Score the theta fit on the FINE partition: a collapsed fine lineage L is
    predicted by its retained ancestor's modeled frequency, split uniformly among
    the co-collapsed fine lineages present that day."""
    ll = 0.0
    for date, counts in fine_by_date.items():
        present = list(counts)
        retained = {L: grouping.get(L, L) for L in present}
        m = {}
        for L in present:
            m[retained[L]] = m.get(retained[L], 0) + 1
        for L, n in counts.items():
            f = modeled.get(retained[L], {}).get(date, EPS) / m[retained[L]]
            ll += n * math.log(max(f, EPS))
    return -ll / N


def thetas_in(ds_dir):
    return sorted(int(r.split("_")[1]) for r in os.listdir(ds_dir)
                  if r.startswith("thresh_")
                  and os.path.exists(os.path.join(ds_dir, r, "mlr_results.json")))


def main():
    metric_rows, peak_rows = [], []
    lineage_sets = [d for d in DATASETS if DATASETS[d]["kind"] == "lineage"]
    for dataset in sorted(lineage_sets):
        ds_dir = os.path.join(SCRATCH, dataset)
        if not os.path.isdir(ds_dir):
            continue
        thetas = thetas_in(ds_dir)
        if not thetas:
            continue
        current = DATASETS[dataset]["current"]

        # Fine reference = lowest-theta fit (finest partition).
        fine_dir = os.path.join(ds_dir, f"thresh_{min(thetas)}")
        fine_by_date, fine_totals, N = read_counts(fine_dir)
        fine_mlr = json.load(open(os.path.join(fine_dir, "mlr_results.json")))
        fine_ga = variant_growth_advantages(fine_mlr)  # missing -> 1.0
        fine_lineages = [v for v in fine_totals if v != "other"]
        dmin, dgen = DATASETS[dataset]["min_date"], DATASETS[dataset]["gen"]

        for theta in thetas:
            run_dir = os.path.join(ds_dir, f"thresh_{theta}")
            mlr = json.load(open(os.path.join(run_dir, "mlr_results.json")))
            modeled = variant_modeled_frequencies(mlr)
            ga = variant_growth_advantages(mlr)
            by_date, totals, _ = read_counts(run_dir)
            se = se_proxies(by_date, dmin, dgen)

            widths = []  # (count, SE) over retained lineages
            for v, n in totals.items():
                if v == "other":
                    continue
                s = se.get(v, float("inf"))
                if not math.isfinite(s):
                    continue
                widths.append((n, s))
                peak_rows.append(dict(dataset=dataset, theta=theta, lineage=v,
                                      count=int(n), freq=n / N,
                                      ga=ga.get(v, 1.0), se_proxy=s))

            # over-collapse from the fine reference: re-group fine lineages at theta
            grouping = collapse_grouping(fine_totals, theta)
            groups = {}
            for L in fine_lineages:
                groups.setdefault(grouping.get(L, L), []).append(L)
            signal_lost = 0.0
            max_fit_err = 0.0
            for members in groups.values():
                w = {L: fine_totals[L] / N for L in members}
                wsum = sum(w.values())
                if wsum <= 0:
                    continue
                gmean = sum(w[L] * fine_ga.get(L, 1.0) for L in members) / wsum
                for L in members:
                    dev = abs(fine_ga.get(L, 1.0) - gmean)
                    signal_lost += w[L] * dev * dev
                    if len(members) > 1:  # only collapsed groups contribute error
                        max_fit_err = max(max_fit_err, dev)

            med = sorted(w for _, w in widths)
            median_se = med[len(med) // 2] if med else float("nan")
            wtot = sum(n for n, _ in widths) or 1
            wmean_se = sum(n * w for n, w in widths) / wtot if widths else float("nan")

            n_variants = len(set(modeled) | set(totals))
            metric_rows.append(dict(
                dataset=dataset, theta=theta, theta_freq=theta / N, is_current=(theta == current),
                N=N, n_variants=n_variants, n_retained_lineages=len(widths),
                coarse_ce=coarse_ce(by_date, totals, N, modeled),
                fixed_ce=fixed_ce(fine_by_date, N, modeled, grouping),
                median_se=median_se, wmean_se=wmean_se,
                signal_lost=signal_lost, max_fit_err=max_fit_err))
            print(f"  {dataset} theta={theta}: n_lin={len(widths)} "
                  f"med_se={median_se:.4f} wmean_se={wmean_se:.4f} "
                  f"signal_lost={signal_lost:.2e} coarse_ce={metric_rows[-1]['coarse_ce']:.3f} "
                  f"fixed_ce={metric_rows[-1]['fixed_ce']:.3f}", flush=True)

    os.makedirs(RESDIR, exist_ok=True)
    mcols = ["dataset", "theta", "theta_freq", "is_current", "N", "n_variants",
             "n_retained_lineages", "coarse_ce", "fixed_ce", "median_se",
             "wmean_se", "signal_lost", "max_fit_err"]
    with open(os.path.join(RESDIR, "lineage_sweep_metrics.tsv"), "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=mcols, delimiter="\t")
        w.writeheader()
        for r in sorted(metric_rows, key=lambda r: (r["dataset"], r["theta"])):
            w.writerow({c: r[c] for c in mcols})
    pcols = ["dataset", "theta", "lineage", "count", "freq", "ga", "se_proxy"]
    with open(os.path.join(RESDIR, "lineage_peaks.tsv"), "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=pcols, delimiter="\t")
        w.writeheader()
        for r in sorted(peak_rows, key=lambda r: (r["dataset"], r["theta"], -r["count"])):
            w.writerow({c: r[c] for c in pcols})
    print(f"\nWrote {len(metric_rows)} metric rows, {len(peak_rows)} lineage rows to {RESDIR}")


if __name__ == "__main__":
    main()
