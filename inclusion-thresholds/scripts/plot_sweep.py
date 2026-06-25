#!/usr/bin/env python3
"""Figures for the clade_min_seq threshold analysis (matplotlib -> figures/*.png).

1. ce_confound.png   -- coarse vs fixed-partition cross-entropy vs threshold:
                        coarse CE *falls* as clades collapse (the metric is
                        confounded by partition coarseness); fixed-partition CE
                        turns up once a real clade is folded into "other".
2. overcollapse.png  -- peak smoothed frequency of the largest collapsed clade
                        vs threshold, as a COUNT and as a FREQUENCY (theta/N).
                        The frequency axis aligns the datasets; the count axis
                        does not. Current production thresholds marked.
3. clade_peaks.png   -- per-clade smoothed peak frequency vs (a) window-total
                        count and (b) mean frequency (count/N). A horizontal
                        f* line is a clean keep/drop rule; no vertical count cut
                        separates real clades from artifacts across datasets.

Prints a recommendation table (current threshold vs the f*-implied threshold).
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESDIR = os.path.join(REPO, "inclusion-thresholds", "results")
FIGDIR = os.path.join(REPO, "inclusion-thresholds", "figures")
COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860"]
FSTAR = 0.01    # peak-frequency boundary between real clades and artifacts (1%)
FMEAN = 0.001   # recommended clade_min_freq: mean-frequency threshold (0.1%)


def datasets(df):
    return list(dict.fromkeys(df["dataset"]))


def fig_ce_confound(df):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for i, ds in enumerate(datasets(df)):
        g = df[df["dataset"] == ds].sort_values("theta")
        c = COLORS[i % len(COLORS)]
        ax1.plot(g["theta"], g["ce_coarse"], "-o", color=c, ms=4, label=ds)
        ax2.plot(g["theta"], g["ce_fixed"], "-o", color=c, ms=4, label=ds)
        cur = g[g["is_current"]]
        if len(cur):
            ax1.plot(cur["theta"], cur["ce_coarse"], "s", color=c, ms=9, mfc="none", mew=1.8)
            ax2.plot(cur["theta"], cur["ce_fixed"], "s", color=c, ms=9, mfc="none", mew=1.8)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.set_xlabel("clade_min_count threshold (count)")
    ax1.set_ylabel("cross-entropy (nats / sequence)")
    ax1.set_title("COARSE CE — scored on each fit's own labels\n(falls as you collapse: confounded)")
    ax2.set_title("FIXED CE — scored on a common label set\n(rises once a real clade is collapsed)")
    ax2.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "ce_confound.png"), dpi=130)
    plt.close(fig)


def fig_overcollapse(df):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for i, ds in enumerate(datasets(df)):
        g = df[df["dataset"] == ds].sort_values("theta")
        c = COLORS[i % len(COLORS)]
        ax1.plot(g["theta"], 100 * g["max_collapsed_peak"], "-o", color=c, ms=4, label=ds)
        ax2.plot(100 * g["theta_freq"], 100 * g["max_collapsed_peak"], "-o", color=c, ms=4, label=ds)
        for ax, x in ((ax1, g["theta"]), (ax2, 100 * g["theta_freq"])):
            cur = g[g["is_current"]]
            if len(cur):
                xc = cur["theta"] if ax is ax1 else 100 * cur["theta_freq"]
                ax.plot(xc, 100 * cur["max_collapsed_peak"], "s", color=c, ms=9, mfc="none", mew=1.8)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.axhline(100 * FSTAR, color="grey", ls="--", lw=1)
        ax.set_ylabel("peak freq. of largest collapsed clade (%)")
    ax1.set_xlabel("threshold as a COUNT")
    ax2.set_xlabel("threshold as a FREQUENCY  theta / N  (%)")
    ax1.set_title("Over-collapse vs count\n(datasets do not align; □ = current)")
    ax2.set_title("Over-collapse vs frequency\n(datasets align; dashed = f* = 1%)")
    ax1.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "overcollapse.png"), dpi=130)
    plt.close(fig)


def fig_clade_peaks(cp):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for i, ds in enumerate(datasets(cp)):
        g = cp[cp["dataset"] == ds]
        c = COLORS[i % len(COLORS)]
        ax1.scatter(g["window_total"], 100 * g["smoothed_peak"], s=22, color=c, label=ds)
        ax2.scatter(100 * g["mean_freq"], 100 * g["smoothed_peak"], s=22, color=c, label=ds)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.axhline(100 * FSTAR, color="grey", ls="--", lw=1)
        ax.set_ylabel("clade smoothed peak frequency (%)")
    ax1.set_xlabel("clade window-total count")
    ax1.set_title("Peak vs COUNT\n(no vertical cut separates real from artifact)")
    ax2.axvline(100 * FMEAN, color="black", ls=":", lw=1.2)
    ax2.set_xlabel("clade mean frequency  count / N  (%)")
    ax2.set_title("Peak vs mean FREQUENCY\n(dashed f*=1% peak; dotted clade_min_freq=0.1%)")
    ax1.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "clade_peaks.png"), dpi=130)
    plt.close(fig)


def recommendation(cp):
    """Per dataset: the count threshold implied by 'keep iff smoothed peak >= f*'.

    Reported as the interval [max window-total among dropped clades (peak<f*),
    min window-total among kept clades (peak>=f*)]; when the upper end is below
    the lower the two classes are cleanly count-separable, otherwise count is an
    ambiguous proxy (an inversion -- a transient real clade with fewer sequences
    than a persistent artifact).
    """
    print(f"\nKeep-rule: smoothed peak >= {100*FSTAR:.1f}%")
    print(f"{'dataset':24} {'N':>9} {'keep_min_ct':>11} {'drop_max_ct':>11} "
          f"{'keep_min_freq%':>14} {'inversion?':>10}")
    for ds in datasets(cp):
        g = cp[cp["dataset"] == ds]
        keep = g[g["smoothed_peak"] >= FSTAR]
        drop = g[g["smoothed_peak"] < FSTAR]
        keep_min = int(keep["window_total"].min()) if len(keep) else 0
        drop_max = int(drop["window_total"].max()) if len(drop) else 0
        N = int(g["N"].iloc[0])
        inv = "YES" if drop_max >= keep_min else "no"
        print(f"{ds:24} {N:>9,} {keep_min:>11} {drop_max:>11} "
              f"{100*keep_min/N:>14.3f} {inv:>10}")


def selection(cp):
    """Choose clade_min_freq: the mean-frequency cut that best reproduces the
    peak>=f* keep/drop rule across all clades, and its per-dataset effect."""
    print(f"\nclade_min_freq selection — misclassification vs the peak>={100*FSTAR:.0f}% rule:")
    print(f"{'cut (%)':>8} {'real_dropped':>13} {'artifact_kept':>14} {'total':>7}")
    for f in (0.0005, 0.0008, 0.001, 0.0015, 0.002):
        rd = ((cp["smoothed_peak"] >= FSTAR) & (cp["mean_freq"] < f)).sum()
        ak = ((cp["smoothed_peak"] < FSTAR) & (cp["mean_freq"] >= f)).sum()
        mark = "  <- clade_min_freq" if abs(f - FMEAN) < 1e-9 else ""
        print(f"{100*f:>8.2f} {rd:>13} {ak:>14} {rd+ak:>7}{mark}")
    print(f"\nclade_min_freq = {FMEAN} ({100*FMEAN:.1f}%) effect per dataset:")
    print(f"{'dataset':24} {'N':>9} {'count=f*N':>10} {'current':>8} {'kept':>5} "
          f"{'max_collapsed_peak%':>19}")
    for ds in datasets(cp):
        g = cp[cp["dataset"] == ds]
        N = int(g["N"].iloc[0])
        thr = FMEAN * N
        dropped = g[g["window_total"] < thr]
        mcp = 100 * dropped["smoothed_peak"].max() if len(dropped) else 0.0
        print(f"{ds:24} {N:>9,} {thr:>10.0f} {'-':>8} {(g['window_total']>=thr).sum():>5} "
              f"{mcp:>19.2f}")


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    df = pd.read_csv(os.path.join(RESDIR, "sweep_metrics.tsv"), sep="\t").sort_values(["dataset", "theta"])
    cp = pd.read_csv(os.path.join(RESDIR, "clade_peaks.tsv"), sep="\t")
    fig_ce_confound(df)
    fig_overcollapse(df)
    fig_clade_peaks(cp)
    recommendation(cp)
    selection(cp)
    print(f"\nFigures -> {FIGDIR}")


if __name__ == "__main__":
    main()
