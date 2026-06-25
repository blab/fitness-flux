#!/usr/bin/env python3
"""Score each swept MLR fit by the cross-entropy of its modeled clade
frequencies against the observed daily counts, plus a complexity penalty (BIC).

Two cross-entropies are reported, and the contrast between them is the point of
the analysis:

* COARSE CE -- the metric originally proposed: score the fit on its OWN label
  set (surviving clades + "other"), i.e. CE = -(1/N) sum n_{v,d} log f_{v,d}
  over the variants that survive threshold theta. This is confounded: raising
  theta merges rare clades into a high-frequency "other" bucket, and a coarser
  partition has lower entropy, so COARSE CE mechanically *improves* as you
  collapse. It therefore has no honest interior optimum (for H3N2 it is
  minimised by collapsing everything to pivot + "other").

* FIXED CE -- score every fit on a COMMON, fine label set (the clades that
  survive the lowest threshold for that dataset). A clade that is collapsed at
  theta is predicted through "other", whose modeled mass is split uniformly
  among the collapsed clades present that day (max-entropy: the coarse model
  genuinely cannot tell them apart). Now collapsing a *real* clade is penalised,
  so FIXED CE has the expected shape and FIXED BIC = 2*N*CE_fixed + k*log(N),
  k = 2*(n_variants-1), has an interior minimum -- the principled threshold.

Also reports `max_collapsed_peak`: the peak *smoothed* frequency (the fit's
`weekly_raw_freq`, i.e. the raw-freq-window empirical average -- NOT raw daily
counts, whose peak is 1.0 on any singleton day) of the single largest clade
folded into "other" at theta. This is the direct over-collapse signal, is purely
empirical (no model needed), and is the basis for the recommended peak-frequency
threshold.

Writes inclusion-thresholds/results/sweep_metrics.tsv (per-threshold metrics)
and inclusion-thresholds/results/clade_peaks.tsv (per-clade window-total count,
mean frequency, and smoothed peak frequency -- the raw material for the
peak-frequency recommendation and the count-vs-frequency contrast).
"""

import csv
import json
import math
import os
import sys

import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SCRATCH = os.path.join(REPO, "inclusion-thresholds", "scratch")
RESULTS = os.path.join(REPO, "inclusion-thresholds", "results", "sweep_metrics.tsv")
sys.path.insert(0, os.path.join(REPO, "fitness-flux-analysis", "scripts"))
sys.path.insert(0, os.path.dirname(__file__))

from ff_io import variant_modeled_frequencies, variant_weekly_frequencies  # noqa: E402
from sweep_clade_min_seq import DATASETS  # noqa: E402

EPS = 1e-6


def fine_reference(ds_dir, thetas):
    """Observed counts and smoothed peaks on the finest label set (lowest theta).

    Returns (by_date, Nd, N, window_total, smoothpeak):
      by_date[date]   = {clade: count}   (the common, theta-invariant outcome space)
      Nd[date]        = daily total
      N               = grand total over the window
      window_total[c] = clade's total count over the window (== the clade_min_seq
                        lever: prepare-data collapses c iff window_total[c] < theta)
      smoothpeak[c]   = peak of the fit's weekly_raw_freq for c (raw-freq-window
                        smoothed empirical frequency; robust unlike raw daily)
    """
    run = os.path.join(ds_dir, f"thresh_{min(thetas)}")
    df = pd.read_csv(os.path.join(run, "prepared_seq_counts.tsv"), sep="\t")
    df["date"] = pd.to_datetime(df["date"]).dt.strftime("%Y-%m-%d")
    df = df.groupby(["date", "variant"], as_index=False)["sequences"].sum()
    by_date, Nd = {}, {}
    for date, sub in df.groupby("date"):
        counts = dict(zip(sub["variant"], sub["sequences"]))
        by_date[date] = counts
        Nd[date] = sum(counts.values())
    N = int(df["sequences"].sum())
    window_total = df.groupby("variant")["sequences"].sum().astype(int).to_dict()

    with open(os.path.join(run, "mlr_results.json")) as h:
        weekly = variant_weekly_frequencies(json.load(h))
    smoothpeak = {c: max(s.values()) for c, s in weekly.items() if s}
    return by_date, Nd, N, window_total, smoothpeak


def score_run(run_dir, by_date, Nd, N, smoothpeak):
    with open(os.path.join(run_dir, "mlr_results.json")) as h:
        mlr = json.load(h)
    modeled = variant_modeled_frequencies(mlr)
    surviving = set(modeled) - {"other"}
    other = modeled.get("other", {})

    ll_coarse = 0.0
    ll_fixed = 0.0
    for date, counts in by_date.items():
        f_other = max(other.get(date, EPS), EPS)
        collapsed_present = [c for c in counts if c not in surviving]
        m = len(collapsed_present) or 1
        # coarse "other" count = sequences of all collapsed clades that day
        n_other = sum(counts[c] for c in collapsed_present)
        if n_other:
            ll_coarse += n_other * math.log(f_other)
        for clade, n in counts.items():
            if clade in surviving:
                f = max(modeled.get(clade, {}).get(date, EPS), EPS)
                ll_coarse += n * math.log(f)
                ll_fixed += n * math.log(f)
            else:
                ll_fixed += n * math.log(max(f_other / m, EPS))
    ce_coarse = -ll_coarse / N
    ce_fixed = -ll_fixed / N

    n_variants = len(surviving) + 1  # + "other"
    k = 2 * (n_variants - 1)
    bic_coarse = 2 * N * ce_coarse + k * math.log(N)
    bic_fixed = 2 * N * ce_fixed + k * math.log(N)

    collapsed = [c for c in smoothpeak if c != "other" and c not in surviving]
    max_collapsed_peak = max((smoothpeak[c] for c in collapsed), default=0.0)
    other_peak = max(other.values()) if other else 0.0
    return dict(N=N, n_variants=n_variants, k=k,
                ce_coarse=ce_coarse, bic_coarse=bic_coarse,
                ce_fixed=ce_fixed, bic_fixed=bic_fixed,
                max_collapsed_peak=max_collapsed_peak, other_peak_freq=other_peak)


def thetas_in(ds_dir):
    return sorted(int(r.split("_")[1]) for r in os.listdir(ds_dir)
                  if r.startswith("thresh_")
                  and os.path.exists(os.path.join(ds_dir, r, "mlr_results.json")))


def main():
    rows = []
    clade_rows = []
    for dataset in sorted(os.listdir(SCRATCH)):
        ds_dir = os.path.join(SCRATCH, dataset)
        if not os.path.isdir(ds_dir):
            continue
        thetas = thetas_in(ds_dir)
        if not thetas:
            continue
        current = DATASETS.get(dataset, {}).get("current")
        by_date, Nd, N, window_total, smoothpeak = fine_reference(ds_dir, thetas)
        for clade, total in sorted(window_total.items()):
            if clade == "other":
                continue
            clade_rows.append(dict(dataset=dataset, clade=clade, N=N,
                                   window_total=total, mean_freq=total / N,
                                   smoothed_peak=smoothpeak.get(clade, 0.0)))
        for theta in thetas:
            run_dir = os.path.join(ds_dir, f"thresh_{theta}")
            try:
                m = score_run(run_dir, by_date, Nd, N, smoothpeak)
            except Exception as e:  # noqa: BLE001
                print(f"  ! {dataset} theta={theta}: {e}", flush=True)
                continue
            m.update(dataset=dataset, theta=theta,
                     theta_freq=theta / m["N"], is_current=(theta == current))
            rows.append(m)
            print(f"  {dataset} theta={theta}: CEc={m['ce_coarse']:.4f} "
                  f"CEf={m['ce_fixed']:.4f} BICf={m['bic_fixed']:.0f} "
                  f"n_var={m['n_variants']} maxcollapse={m['max_collapsed_peak']:.3f}",
                  flush=True)

    cols = ["dataset", "theta", "theta_freq", "is_current", "N", "n_variants", "k",
            "ce_coarse", "bic_coarse", "ce_fixed", "bic_fixed",
            "max_collapsed_peak", "other_peak_freq"]
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["dataset"], r["theta"])):
            w.writerow({c: r[c] for c in cols})
    print(f"\nWrote {len(rows)} rows to {RESULTS}")

    peaks_path = os.path.join(os.path.dirname(RESULTS), "clade_peaks.tsv")
    ccols = ["dataset", "clade", "N", "window_total", "mean_freq", "smoothed_peak"]
    with open(peaks_path, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=ccols, delimiter="\t")
        w.writeheader()
        for r in sorted(clade_rows, key=lambda r: (r["dataset"], -r["window_total"])):
            w.writerow({c: r[c] for c in ccols})
    print(f"Wrote {len(clade_rows)} clade rows to {peaks_path}")


if __name__ == "__main__":
    main()
