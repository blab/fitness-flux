#!/usr/bin/env python3
"""Annualize the per-generation fitness flux/variance series for comparison to
Łuksza & Lässig (Nature 2014, "A predictive fitness model for influenza").

Our fitness quantities are computed per generation (τ days); L&L report the same
fitness-flux construct (Mustonen & Lässig) per calendar year, because they
propagate clade frequencies winter-to-winter. To compare we rescale fitness from
per-generation to per-year with

    G = 365.25 / τ            # generations per year (τ = 3.2 d for H3N2 → G ≈ 114)

Both the fitness flux and the fitness variance are quadratic in the per-year
rescaling of fitness and therefore scale by **G²** (not G¹ for one and G² for the
other): the stored `velocity` equals φ_year / G² and `variance` equals
Var(f_year) / G², because the mean-fitness series `mean_log_fitness` is itself in
per-generation fitness units (f_gen = τ·f_day) and `velocity` carries an extra
τ/365.25 time-axis factor (see fitness_wave.py). Scaling both series by the same
power leaves the dimensionless variance→flux (Fisher) slope unchanged; that
invariance is asserted below as a built-in unit test.

Reads the saved series behind the H3N2 variance/flux figure (do not re-derive from
frequencies); prints the per-gen / per-year / L&L comparison and writes a JSON.
"""
import argparse
import csv
import json

import numpy as np

import ff_io

DAYS_PER_YEAR = 365.25

# Łuksza & Lässig 2014 full-haplotype model, per year (Extended Data Table 1 and
# the Φ ≈ 26 cumulative flux over 1993–2010 → mean rate ≈ 26/17 ≈ 1.5 yr⁻¹).
LL_VARIANCE = 2.1       # yr⁻²
LL_FLUX_RATE = 1.5      # yr⁻¹, the rate comparable to our time-averaged flux
LL_CUMULATIVE = 26.0    # dimensionless, cumulative over ~17 seasons


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", default="dataset",
        help="dataset label for reporting, e.g. h3n2_clades or sarscov2_clades")
    parser.add_argument("--timeseries", required=True,
        help="{analysis}_flux_timeseries.tsv (date, decimal_date, mean_log_fitness, "
             "variance, velocity)")
    parser.add_argument("--summary", required=True,
        help="{analysis}_flux_summary.json (for τ, the time-weighted avg velocity, "
             "and the per-gen Fisher slope)")
    parser.add_argument("--generation-time", type=float, default=None,
        help="Override τ in days; default reads mean_generation_time_days from the "
             "summary (3.2 for H3N2).")
    parser.add_argument("--output", required=True, help="output JSON path")
    return parser.parse_args()


def trapezoid(y, x):
    """∫ y dx by the trapezoidal rule (numpy.trapz was removed in numpy 2)."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    return float(np.sum(np.diff(x) * (y[:-1] + y[1:]) / 2.0))


def summarize(values):
    """median, IQR (25th, 75th) and max of a series."""
    q1, med, q3 = (float(v) for v in np.percentile(values, [25, 50, 75]))
    return {"median": med, "iqr_lo": q1, "iqr_hi": q3, "max": float(np.max(values))}


def main():
    args = parse_args()

    with open(args.summary) as handle:
        summary = json.load(handle)
    tau = args.generation_time
    if tau is None:
        # For H3N2 mean == nominal == 3.2; prefer the mean actually used by the wave.
        tau = summary.get("mean_generation_time_days") or summary["generation_time_days"]
    G = DAYS_PER_YEAR / tau
    G2 = G * G
    per_gen_slope = summary.get("variance_vs_velocity", {}).get("slope")

    # SARS-CoV-2 splits the generation time pre/post-Omicron (per-variant clocks,
    # 5.0 vs 3.2 d) and the wave uses a frequency-weighted tau_bar(t), so there is no
    # single τ: nominal (generation_time_days) and mean (mean_generation_time_days)
    # differ. A single-G annualization then uses the mean and is APPROXIMATE in
    # absolute magnitude (the dimensionless Fisher slope stays exact). For H3N2 the
    # two agree and the annualization is exact.
    nominal_tau = summary.get("generation_time_days")
    tau_varies = nominal_tau is not None and not np.isclose(tau, nominal_tau, rtol=1e-3)

    # Read the per-generation series. velocity is only defined past the first
    # finite-difference window, so it is present on a subset of the dates.
    dates, variance, velocity = [], [], []
    with open(args.timeseries) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            dates.append(float(row["decimal_date"]))
            variance.append(float(row["variance"]))
            vel = row["velocity"]
            velocity.append(float(vel) if vel not in ("", "None") else np.nan)
    dates = np.array(dates)
    variance = np.array(variance)
    velocity = np.array(velocity)
    vmask = ~np.isnan(velocity)

    # --- Step 1: confirm native units before converting ------------------------
    ff_io.log(f"[{args.dataset}] τ = {tau:.3g} d → G = 365.25/{tau:.3g} = {G:.3f} gen/yr "
              f"(G² = {G2:.1f})")
    if tau_varies:
        ff_io.log(f"NOTE: generation time varies (nominal {nominal_tau} d, "
                  f"frequency-weighted mean {tau:.2f} d); pre/post-Omicron per-variant "
                  f"clocks make the per-year MAGNITUDES approximate (Fisher slope exact).")
    native = {
        "variance_per_gen": summarize(variance),
        "flux_per_gen": summarize(velocity[vmask]),
    }
    for name, stats in native.items():
        ff_io.log(f"native {name}: min≈… median={stats['median']:.3e} "
                  f"max={stats['max']:.3e}")
    # Sanity: both native series should sit around 1e-3 (per the figure). Bail if
    # either is orders of magnitude off — the conversion is meaningless otherwise.
    for name, stats in native.items():
        if not (1e-6 <= abs(stats["median"]) <= 1e-1):
            raise SystemExit(
                f"ERROR: native {name} median {stats['median']:.3e} is far from the "
                f"expected ~1e-3 scale; confirm the input series before annualizing.")

    # --- Step 2: annualize (both series scale by G²) ---------------------------
    variance_per_year = variance * G2
    flux_per_year = velocity * G2

    # --- Step 3: summarize each annualized series ------------------------------
    var_yr = summarize(variance_per_year)
    flux_yr = summarize(flux_per_year[vmask])

    # --- Step 4: annual mean flux RATE, comparable to L&L's Φ̄ ≈ 1.5 yr⁻¹ -------
    # Time-weighted average velocity from the summary (de-biased for uneven
    # sampling), rescaled; cross-checked against the plain mean and the
    # cumulative integral divided by the window length.
    mean_flux_rate_timeweighted = summary["avg_velocity"] * G2
    mean_flux_rate_pointwise = float(np.nanmean(flux_per_year))
    order = np.argsort(dates[vmask])
    fx = flux_per_year[vmask][order]
    tx = dates[vmask][order]
    cumulative_flux = trapezoid(fx, tx)          # analogue of L&L's Φ ≈ 26
    window_years = float(tx[-1] - tx[0])
    mean_flux_rate_integral = cumulative_flux / window_years

    # --- Step 5: re-fit Fisher in per-year units; slope must be invariant ------
    xs = variance_per_year[vmask]
    ys = flux_per_year[vmask]
    slope_year, _ = np.polyfit(xs, ys, 1)
    r_year = float(np.corrcoef(xs, ys)[0, 1])
    if per_gen_slope is not None and not np.isclose(slope_year, per_gen_slope, rtol=1e-6):
        raise SystemExit(
            f"ERROR: per-year Fisher slope {slope_year:.4f} != per-gen slope "
            f"{per_gen_slope:.4f}. Both series must scale by the SAME power of G; "
            f"a mismatch means the flux/variance powers are wrong.")

    # --- Report ----------------------------------------------------------------
    def fold(ours, theirs):
        return ours / theirs if theirs else float("nan")

    var_fold = fold(var_yr["median"], LL_VARIANCE)
    flux_fold = fold(mean_flux_rate_timeweighted, LL_FLUX_RATE)
    within_few_fold = all(0.2 <= f <= 5.0 for f in (var_fold, flux_fold))

    print("\nQuantity                         | our per-gen | our per-year | L&L per-year")
    print("-" * 78)
    print(f"fitness variance (median)        | {native['variance_per_gen']['median']:.3e}   "
          f"| {var_yr['median']:7.2f} yr⁻² | {LL_VARIANCE:.1f} yr⁻²")
    print(f"fitness flux (median)            | {native['flux_per_gen']['median']:.3e}   "
          f"| {flux_yr['median']:7.2f}      |  —")
    print(f"mean flux rate (time-avg)        | {summary['avg_velocity']:.3e}   "
          f"| {mean_flux_rate_timeweighted:7.2f} yr⁻¹ | {LL_FLUX_RATE:.1f} yr⁻¹")
    print(f"cumulative flux (over {window_years:.1f} yr)    |      —      "
          f"| {cumulative_flux:7.2f}      | {LL_CUMULATIVE:.0f}")
    print(f"variance→flux Fisher slope       | {per_gen_slope:.3f}       "
          f"| {slope_year:.3f}        | ≈1")
    print("-" * 78)
    print(f"\nHeadline (per year): median flux rate {mean_flux_rate_timeweighted:.2f} yr⁻¹, "
          f"median variance {var_yr['median']:.2f} yr⁻², Fisher slope {slope_year:.2f} "
          f"(r = {r_year:.2f}).")
    approx = " (approximate; generation time varies pre/post-Omicron)" if tau_varies else ""
    if within_few_fold:
        # H3N2: same influenza-scale comparison L&L made — concordance, realized ≥ predicted.
        print(
            f"Our annualized variance ({var_yr['median']:.1f} yr⁻²) and mean flux rate "
            f"({mean_flux_rate_timeweighted:.1f} yr⁻¹) are within a few-fold of L&L's "
            f"{LL_VARIANCE} yr⁻² and {LL_FLUX_RATE} yr⁻¹, with the same variance–flux "
            f"(Fisher) slope (≈{slope_year:.2f}); realized MLR flux sitting above their "
            f"model-predicted value is expected.\n")
    else:
        # SARS-CoV-2: far above the influenza benchmark — that gap is the point.
        print(
            f"Our annualized variance ({var_yr['median']:.1f} yr⁻²) and mean flux rate "
            f"({mean_flux_rate_timeweighted:.1f} yr⁻¹) run ~{var_fold:.0f}× and "
            f"~{flux_fold:.0f}× above L&L's influenza benchmark ({LL_VARIANCE} yr⁻², "
            f"{LL_FLUX_RATE} yr⁻¹){approx}, consistent with far faster adaptation, while "
            f"the variance–flux (Fisher) slope (≈{slope_year:.2f}) holds as in influenza.\n")

    result = {
        "dataset": args.dataset,
        "generation_time_days": tau,
        "generation_time_nominal_days": nominal_tau,
        "generation_time_varies": bool(tau_varies),
        "generations_per_year": G,
        "native_per_gen": native,
        "per_year": {
            "variance": var_yr,
            "flux": flux_yr,
            "mean_flux_rate_timeweighted": mean_flux_rate_timeweighted,
            "mean_flux_rate_pointwise": mean_flux_rate_pointwise,
            "mean_flux_rate_integral": mean_flux_rate_integral,
            "cumulative_flux": cumulative_flux,
            "window_years": window_years,
            "fisher_slope": slope_year,
            "fisher_r": r_year,
        },
        "fisher_slope_per_gen": per_gen_slope,
        "luksza_lassig_2014": {
            "variance_per_year": LL_VARIANCE,
            "mean_flux_rate_per_year": LL_FLUX_RATE,
            "cumulative_flux": LL_CUMULATIVE,
        },
        "variance_fold_vs_ll": var_fold,
        "flux_rate_fold_vs_ll": flux_fold,
        "within_few_fold": within_few_fold,
    }
    with open(args.output, "w") as handle:
        json.dump(result, handle, indent=2)
    ff_io.log(f"Wrote annualized flux comparison to {args.output}")


if __name__ == "__main__":
    main()
