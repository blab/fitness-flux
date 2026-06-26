#!/usr/bin/env python3
"""Figures for the lineage `collapse_threshold` analysis -> figures/lineage_*.png.

1. lineage_noise.png      -- per-lineage growth-rate SE (analytic) vs COUNT and vs
                             FREQUENCY (pooled across windows). Hypothesis: vs count
                             the points collapse onto one ~1/sqrt(n) curve (noise is
                             count-governed); vs frequency they do not.
2. lineage_noise_theta.png-- median SE vs threshold per window (noise growth
                             at low theta).
3. lineage_overcollapse.png-- fitness signal lost (freq-weighted variance of fine
                             lineage growth advantages erased by collapsing) vs
                             threshold, count and frequency axes.
4. lineage_confound.png   -- coarse vs fixed-partition CE (the confound, as in clades).
5. lineage_verdict.png    -- noise floor n* (count where SE crosses a target)
                             and over-collapse onset per window vs N: count vs frequency.

Prints a recommendation table.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESDIR = os.path.join(REPO, "inclusion-thresholds", "results")
FIGDIR = os.path.join(REPO, "inclusion-thresholds", "figures")
COLORS = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860"]
WSTAR = 0.01  # target growth-rate SE (1/gen) noise floor (proxy; NUTS calibrates absolute level)


def datasets(df):
    return list(dict.fromkeys(df["dataset"]))


def fig_noise(cp):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for i, ds in enumerate(datasets(cp)):
        g = cp[cp["dataset"] == ds]
        c = COLORS[i % len(COLORS)]
        ax1.scatter(g["count"], g["se_proxy"], s=12, alpha=0.5, color=c, label=ds)
        ax2.scatter(100 * g["freq"], g["se_proxy"], s=12, alpha=0.5, color=c, label=ds)
    # 1/sqrt(n) guide fit on the pooled cloud
    m = cp[(cp["se_proxy"] > 0) & (cp["count"] > 0)]
    C = float(np.median(m["se_proxy"] * np.sqrt(m["count"])))
    xs = np.array([m["count"].min(), m["count"].max()])
    ax1.plot(xs, C / np.sqrt(xs), "k--", lw=1.2, label=f"{C:.2f}/√n")
    for ax in (ax1, ax2):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.axhline(WSTAR, color="grey", ls=":", lw=1)
        ax.set_ylabel("growth-rate SE (1/gen)")
    ax1.set_xlabel("lineage count n")
    ax1.set_title("Noise vs COUNT\n(collapses onto 1/√n across windows)")
    ax2.set_xlabel("lineage frequency n/N (%)")
    ax2.set_title("Noise vs FREQUENCY\n(does not collapse; dotted = target W*)")
    ax1.legend(fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lineage_noise.png"), dpi=130)
    plt.close(fig)


def fig_noise_theta(df):
    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    for i, ds in enumerate(datasets(df)):
        g = df[df["dataset"] == ds].sort_values("theta")
        c = COLORS[i % len(COLORS)]
        ax.plot(g["theta"], g["median_se"], "-o", color=c, ms=4, label=ds)
        cur = g[g["is_current"]]
        if len(cur):
            ax.plot(cur["theta"], cur["median_se"], "s", color=c, ms=9, mfc="none", mew=1.8)
    ax.axhline(WSTAR, color="grey", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("collapse_threshold (count)")
    ax.set_ylabel("median lineage SE (1/gen)")
    ax.set_title("Noise vs threshold (□ = current 1000; dotted = W*)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lineage_noise_theta.png"), dpi=130)
    plt.close(fig)


def fig_overcollapse(df):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for i, ds in enumerate(datasets(df)):
        g = df[df["dataset"] == ds].sort_values("theta")
        c = COLORS[i % len(COLORS)]
        ax1.plot(g["theta"], g["signal_lost"], "-o", color=c, ms=4, label=ds)
        ax2.plot(100 * g["theta_freq"], g["signal_lost"], "-o", color=c, ms=4, label=ds)
        for ax, x in ((ax1, g["theta"]), (ax2, 100 * g["theta_freq"])):
            cur = g[g["is_current"]]
            if len(cur):
                xc = cur["theta"] if ax is ax1 else 100 * cur["theta_freq"]
                ax.plot(xc, cur["signal_lost"], "s", color=c, ms=9, mfc="none", mew=1.8)
    for ax in (ax1, ax2):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_ylim(bottom=1e-5)  # the lowest-theta point is ~0 (fine ref vs itself)
        ax.set_ylabel("fitness signal lost (freq-weighted var)")
    ax1.set_xlabel("collapse_threshold (count)")
    ax1.set_title("Over-collapse vs COUNT (□ = current)")
    ax2.set_xlabel("threshold as frequency θ/N (%)")
    ax2.set_title("Over-collapse vs FREQUENCY")
    ax1.legend(fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lineage_overcollapse.png"), dpi=130)
    plt.close(fig)


def fig_confound(df):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for i, ds in enumerate(datasets(df)):
        g = df[df["dataset"] == ds].sort_values("theta")
        c = COLORS[i % len(COLORS)]
        ax1.plot(g["theta"], g["coarse_ce"], "-o", color=c, ms=4, label=ds)
        ax2.plot(g["theta"], g["fixed_ce"], "-o", color=c, ms=4, label=ds)
    for ax in (ax1, ax2):
        ax.set_xscale("log"); ax.set_xlabel("collapse_threshold (count)")
    ax1.set_ylabel("cross-entropy (nats / sequence)")
    ax1.set_title("COARSE CE — own labels (rewards collapsing)")
    ax2.set_title("FIXED CE — common fine partition")
    ax2.legend(fontsize=7.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lineage_confound.png"), dpi=130)
    plt.close(fig)


def noise_floor(cp, ds):
    """Smallest count whose binned-median HDI width <= WSTAR (the count noise floor)."""
    g = cp[cp["dataset"] == ds].sort_values("count")
    if g.empty:
        return None
    # rolling median over count-sorted widths
    g = g.assign(roll=g["se_proxy"].rolling(7, min_periods=3, center=True).median())
    ok = g[g["roll"] <= WSTAR]
    return int(ok["count"].min()) if len(ok) else None


def recommendation(df, cp):
    print(f"\nNoise floor: smallest count with median SE <= {WSTAR}")
    print(f"{'dataset':26} {'N':>10} {'noise_floor_ct':>14} {'as % of N':>10} "
          f"{'current=1000 %N':>15} {'signal_lost@1000':>17}")
    for ds in datasets(df):
        g = df[df["dataset"] == ds]
        N = int(g["N"].iloc[0])
        nf = noise_floor(cp, ds)
        cur = g[g["is_current"]]
        sl = float(cur["signal_lost"].iloc[0]) if len(cur) else float("nan")
        nf_s = f"{nf}" if nf else "n/a"
        nf_pct = f"{100*nf/N:.3f}" if nf else "n/a"
        print(f"{ds:26} {N:>10,} {nf_s:>14} {nf_pct:>10} {100*1000/N:>15.3f} {sl:>17.2e}")


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    df = pd.read_csv(os.path.join(RESDIR, "lineage_sweep_metrics.tsv"), sep="\t")
    cp = pd.read_csv(os.path.join(RESDIR, "lineage_peaks.tsv"), sep="\t")
    fig_noise(cp)
    fig_noise_theta(df)
    fig_overcollapse(df)
    fig_confound(df)
    recommendation(df, cp)
    print(f"\nFigures -> {FIGDIR}")


if __name__ == "__main__":
    main()
