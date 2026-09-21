#!/usr/bin/env python3
"""Forecast-accuracy analysis: MLR vs a naive model (Abousamra et al. Fig 2B).

Adaptation of the retrospective forecast evaluation in Abousamra, Figgins &
Bedford 2024 ("Fitness models provide accurate short-term forecasts of
SARS-CoV-2 variant frequency") to this repo's window cadence. It measures the
mean absolute error (MAE) of clade-frequency predictions as a function of
forecast lead time, comparing the MLR model to a naive (persistence) model. It
runs entirely on the already-computed per-window MLR fits (mlr_results.json) and
does not re-fit anything.

Window pairing
--------------
Each dataset is fit in overlapping windows that slide forward by a fixed amount
(SARS-CoV-2: 1-year windows sliding 6 months; flu: 2-year windows sliding 1
year). A window W's successor W+1 is the one fit ~one slide later, found by a
day-gap tolerance on the actual final data dates (--slide-days; NOT relativedelta
arithmetic, which mis-lands on end-of-month dates, and NOT the nominal config
bounds, since sparse flu data can end early). A window whose successor is missing
(e.g. across a bridge window) is skipped.

The forecast horizon is INDEPENDENT of the slide. We treat W's final date T as
the "date of estimation" and evaluate a fixed hindcast (in-window, lead <= 0) and
forecast (lead > 0) window against W+1's empirical smoothed frequency as truth.
For flu the successor spans a full year beyond W, but we still score only
--forecast-days (default 180): beyond ~6 months flu data is sparse and MLR and
naive converge, so the 6-12 month tail is uninformative and excluded.

Predictions (over W's *named* clades S = metadata["variants"] minus "other")
--------------------------------------------------------------------------------
The "other" bucket is excluded and the named clades are renormalized to sum 1, so
we score how well named-clade relative frequencies are predicted (neither model can
forecast the growing "other" bucket of newly emerging clades).
delta_v (per-day log-growth slope) is recovered exactly from the exported growth
advantage: the MLR export writes ga_v = exp(delta_v * tau_v) with the same
per-variant tau_v used at fit time (scripts/generation_time.py), so
delta_v = ln(ga_v) / tau_v. Since evofr's freq_v(t) = softmax_v(intercept + slope*t)
and softmax is shift-invariant, softmax_v[ln freq_v(T) + delta_v*lead]
reconstructs the model's own forward projection.

  MLR:   lead <= 0  -> W's in-window modeled freq_v(t) (site "freq", median)
         lead >  0  -> softmax_v[ ln max(freq_v(T), eps) + delta_v*lead ]
  Naive: lead <= 0  -> identical to MLR (both are the in-window fit)
         lead >  0  -> freq_v(T) held constant (flat forward from T)

The models coincide at lead 0 and diverge only for lead > 0, so the comparison
isolates the forward projection.

Truth and error
---------------
Truth at date t is W+1's weekly_raw_freq restricted to S and renormalized over S
(W+1's "other" bucket and post-W emergent clades are dropped). Absolute error
(Abousamra Eq 2) is AE_t = (1/n) * sum_{v in S} |truth_v(t) - pred_v(t)|, n = |S|.
We aggregate AE across window pairs on a shared lead grid (weekly bins) to get MAE(lead).

Limitations
-----------
- ga / freq are exported rounded to 3 decimals; a clade at freq == 0 at T can
  only grow from the eps floor, so a clade poised to sweep but rounded to 0 is
  under-forecast by both models (inherent to using the median freq export).
- Over a 6-month horizon clades that emerge after W cannot be forecast (they
  land in truth's "other"), so absolute magnitudes are larger than, and not
  directly comparable to, the paper's 30-day numbers; the MLR < naive ordering
  and rough scale are the checks.
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from datetime import date

import ff_io

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from generation_time import generation_time_for, load_aliasor  # noqa: E402

# A window's successor sits ~one slide later (SARS-CoV-2 slides 6 months; flu
# slides 1 year). We accept a successor whose final data date is within +/-20% of
# the slide, which spans the real spread (e.g. 181-184 d for a 6-month slide,
# 324-407 d for a 1-year slide) while rejecting a doubled jump when the true
# successor is absent, and same-year duplicates from bridge windows.
GAP_TOLERANCE = 0.20


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="e.g. sarscov2_clades")
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument(
        "--slide-days",
        type=int,
        default=180,
        help="days ahead to find each window's successor (the truth source); 180 pairs "
        "the window two quarters ahead for a 6-month horizon. Independent of --forecast-days.",
    )
    parser.add_argument(
        "--generation-time",
        type=float,
        default=3.2,
        help="generation time in days (tau; used to recover delta = ln(ga)/tau)",
    )
    parser.add_argument(
        "--generation-time-pre-omicron",
        type=float,
        default=None,
        help="if set, pre-Omicron variants use this tau and others --generation-time",
    )
    parser.add_argument(
        "--variant-classification",
        choices=["clades", "lineages"],
        default=None,
        help="how to read variant names for the pre/post-Omicron split",
    )
    parser.add_argument(
        "--aliasing",
        type=str,
        default=None,
        help="optional local Pango alias_key.json for lineage classification",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-3,
        help="floor for freq/ga in log space (matches the 3-decimal export granularity)",
    )
    parser.add_argument("--lead-bin-days", type=int, default=7)
    parser.add_argument("--hindcast-days", type=int, default=90)
    parser.add_argument("--forecast-days", type=int, default=180)
    parser.add_argument("--detail-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def day_gap(a_iso, b_iso):
    """a - b in days, from ISO YYYY-MM-DD strings."""
    return (date.fromisoformat(a_iso) - date.fromisoformat(b_iso)).days


def load_windows(mlr_dir, dataset):
    """(timepoint, mlr, final_date) per window, sorted by final date."""
    windows = []
    for tp in ff_io.seasonal_timepoints(mlr_dir, dataset):
        mlr = ff_io.load_mlr(mlr_dir, dataset, tp)
        windows.append((tp, mlr, mlr["metadata"]["dates"][-1]))
    windows.sort(key=lambda w: w[2])
    return windows


def build_pairs(windows, slide_days):
    """(W, successor) pairs where the successor's final data date is ~slide_days later.

    Windows are sorted by final date. The successor is not necessarily the adjacent
    window: with a 3-month slide but a 6-month (180 d) forecast horizon, each window
    pairs with the one two quarters ahead. Take the nearest window whose gap falls in
    the tolerance band; a window with no such successor (e.g. across a data gap) is skipped.
    """
    lo = slide_days * (1 - GAP_TOLERANCE)
    hi = slide_days * (1 + GAP_TOLERANCE)
    pairs = []
    for i, cur in enumerate(windows):
        for nxt in windows[i + 1:]:
            gap = day_gap(nxt[2], cur[2])
            if gap > hi:
                break  # sorted by final date: no later window qualifies
            if gap >= lo:
                pairs.append((cur, nxt))
                break  # nearest qualifying successor
    return pairs


def make_tau_fn(args):
    """variant -> generation time tau. Uniform when no pre-Omicron split."""
    if args.generation_time_pre_omicron is None or args.variant_classification is None:
        return lambda v: args.generation_time
    aliasor = (
        load_aliasor(args.aliasing)
        if args.variant_classification == "lineages"
        else None
    )
    return lambda v: generation_time_for(
        v,
        args.variant_classification,
        args.generation_time_pre_omicron,
        args.generation_time,
        aliasor,
    )


def softmax(logits):
    """logits: {variant: value} -> {variant: probability}."""
    hi = max(logits.values())
    exps = {v: math.exp(x - hi) for v, x in logits.items()}
    total = sum(exps.values())
    return {v: e / total for v, e in exps.items()}


def normalize(values):
    total = sum(values.values())
    if total <= 0:
        return {v: 0.0 for v in values}
    return {v: x / total for v, x in values.items()}


def remap_truth(truth_next, variant_set, d):
    """W+1 weekly_raw_freq at date d, restricted to W's named clades (variant_set)
    and renormalized to sum 1 over them.

    Mass outside variant_set — W+1's "other" bucket and clades that emerged after W
    — is dropped, matching the named-clade-only convention (we score how well the
    named clades' relative frequencies are predicted, since neither model can
    forecast the growing "other"). Returns None when d has no non-null named truth.
    """
    s = set(variant_set)
    kept = {}
    for variant, series in truth_next.items():
        if variant not in s:
            continue
        value = series.get(d)
        if value is not None:
            kept[variant] = value
    total = sum(kept.values())
    if not kept or total <= 0:
        return None
    return {v: kept.get(v, 0.0) / total for v in variant_set}


def forecast_pair(cur, nxt, tau_fn, epsilon, hindcast_days, forecast_days, location=None):
    """Rows (date, lead, n, model, abs_error) for one (W, W+1) pair in one region.

    Frequencies are over W's *named* clades (the "other" bucket is excluded and the
    named clades are renormalized to sum 1), for both predictions and truth. All
    series are read for ``location`` (a region); None => primary (single-region).
    """
    _, m, T = cur
    _, nm, _ = nxt
    variant_set = [v for v in m["metadata"]["variants"] if v != "other"]
    n = len(variant_set)
    if n < 2:
        return []  # need >=2 named clades to define relative frequencies

    ga = ff_io.variant_growth_advantages(m, location=location)          # missing / "other" -> 1.0
    modeled = ff_io.variant_modeled_frequencies(m, location=location)   # {v: {date: freq}} (median)
    truth_next = ff_io.variant_weekly_frequencies(nm, location=location)

    freqT = normalize({v: (modeled.get(v, {}).get(T) or 0.0) for v in variant_set})
    delta = {v: math.log(max(ga.get(v, 1.0), epsilon)) / tau_fn(v) for v in variant_set}
    logf0 = {v: math.log(max(freqT[v], epsilon)) for v in variant_set}

    rows = []
    for d in nm["metadata"]["dates"]:
        lead = day_gap(d, T)
        if lead < -hindcast_days or lead > forecast_days:
            continue
        truth = remap_truth(truth_next, variant_set, d)
        if truth is None:
            continue
        if lead <= 0:
            # Hindcast: both models are the in-window fit, so they coincide.
            pred_mlr = normalize({v: (modeled.get(v, {}).get(d) or 0.0) for v in variant_set})
            pred_naive = pred_mlr
        else:
            pred_mlr = softmax({v: logf0[v] + delta[v] * lead for v in variant_set})
            pred_naive = freqT
        ae_mlr = sum(abs(truth[v] - pred_mlr[v]) for v in variant_set) / n
        ae_naive = sum(abs(truth[v] - pred_naive[v]) for v in variant_set) / n
        rows.append((d, lead, n, "mlr", ae_mlr))
        rows.append((d, lead, n, "naive", ae_naive))
    return rows


def aggregate(all_rows, bin_days):
    """MAE(lead) curve: mean abs_error per model within weekly lead bins (pooled over regions)."""
    bucket = defaultdict(lambda: defaultdict(list))  # bin_index -> model -> [ae]
    for _region, _pair, _win, _nxt, _d, lead, _n, model, ae in all_rows:
        bucket[lead // bin_days][model].append(ae)
    curve = []
    for b in sorted(bucket):
        entry = {
            "lead": b * bin_days + (bin_days - 1) / 2.0,
            "lead_lo": b * bin_days,
            "lead_hi": b * bin_days + bin_days - 1,
        }
        for model in ("mlr", "naive"):
            values = bucket[b].get(model, [])
            entry[model] = sum(values) / len(values) if values else None
            entry["n_" + model] = len(values)
        curve.append(entry)
    return curve


def horizon_means(all_rows):
    """Per-model mean AE over the hindcast (lead <= 0) and forecast (lead > 0)."""
    acc = {"hindcast": defaultdict(list), "forecast": defaultdict(list)}
    for _region, _pair, _win, _nxt, _d, lead, _n, model, ae in all_rows:
        acc["forecast" if lead > 0 else "hindcast"][model].append(ae)
    out = {}
    for horizon, by_model in acc.items():
        out[horizon] = {
            model: (sum(v) / len(v) if v else None)
            for model, v in {"mlr": by_model.get("mlr", []), "naive": by_model.get("naive", [])}.items()
        }
    return out


def main():
    args = parse_args()
    tau_fn = make_tau_fn(args)

    windows = load_windows(args.mlr_dir, args.dataset)
    # Regions present across the dataset's windows (each window's hierarchical fit
    # carries one series per sufficiently-sampled region). Forecasts are paired and
    # scored WITHIN a region (W and W+1 from the same region).
    all_regions = sorted(set().union(*[set(ff_io.regions(w[1])) for w in windows])) if windows else []

    all_rows = []  # (region, pair_label, window, next_window, date, lead, n, model, ae)
    pairs_by_region = {}
    for region in all_regions:
        region_windows = [w for w in windows if region in ff_io.regions(w[1])]
        pairs = build_pairs(region_windows, args.slide_days)
        pairs_by_region[region] = [f"{cur[0]}->{nxt[0]}" for cur, nxt in pairs]
        for cur, nxt in pairs:
            label = f"{cur[0]}->{nxt[0]}"
            for d, lead, n, model, ae in forecast_pair(
                cur, nxt, tau_fn, args.epsilon, args.hindcast_days, args.forecast_days,
                location=region,
            ):
                all_rows.append((region, label, cur[0], nxt[0], d, lead, n, model, ae))
        ff_io.log(f"  {region}: {len(region_windows)} windows, {len(pairs)} pairs")

    total_pairs = sum(len(p) for p in pairs_by_region.values())
    ff_io.log(
        f"{args.dataset}: {len(all_regions)} regions, {total_pairs} region-window pairs "
        f"(slide {args.slide_days}d, horizon +{args.forecast_days}d / -{args.hindcast_days}d)"
    )

    with open(args.detail_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["region", "window", "next_window", "date", "lead_days", "n_variants", "model", "abs_error"]
        )
        for region, _label, window, next_window, d, lead, n, model, ae in all_rows:
            writer.writerow([region, window, next_window, d, lead, n, model, f"{ae:.6f}"])
    ff_io.log(f"Wrote {len(all_rows)} detail rows to {args.detail_output}")

    by_region = {
        region: {"n_pairs": len(pairs_by_region[region]),
                 **horizon_means([r for r in all_rows if r[0] == region])}
        for region in all_regions
    }

    summary = {
        "dataset": args.dataset,
        "epsilon": args.epsilon,
        "lead_bin_days": args.lead_bin_days,
        "hindcast_days": args.hindcast_days,
        "forecast_days": args.forecast_days,
        "regions": all_regions,
        "n_pairs": total_pairs,
        "pairs_by_region": pairs_by_region,
        "reference": 0.05,
        "overall": horizon_means(all_rows),
        "by_region": by_region,
        "curve": aggregate(all_rows, args.lead_bin_days),
    }
    with open(args.summary_output, "w") as handle:
        json.dump(summary, handle, indent=2)
    ff_io.log(f"Wrote forecast-accuracy summary to {args.summary_output}")


if __name__ == "__main__":
    main()
