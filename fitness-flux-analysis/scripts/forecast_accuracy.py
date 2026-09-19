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
Each dataset (e.g. sarscov2_clades) is fit in overlapping windows that slide
forward ~6 months, so a window W's successor W+1 is fit on data shifted +6
months and extends exactly 6 months beyond W. We treat W's final date T as the
"date of estimation" and evaluate a hindcast (in-window, lead <= 0) plus a
forecast (lead > 0) against W+1's empirical smoothed frequency as retrospective
truth. Successors are found by a day-gap tolerance on the windows' final dates
(NOT relativedelta arithmetic, which mis-lands on end-of-month dates); a window
whose +6-month successor is missing from config["datasets"] is skipped.

Predictions (over W's fitted variant set S = metadata["variants"], incl. "other")
--------------------------------------------------------------------------------
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
Truth at date t is W+1's weekly_raw_freq remapped onto S (any W+1 variant not in
S folded into "other"), renormalized over S. Absolute error (Abousamra Eq 2) is
AE_t = (1/n) * sum_{v in S} |truth_v(t) - pred_v(t)|, n = |S|. We aggregate AE
across window pairs on a shared lead grid (weekly bins) to get MAE(lead).

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

# Successor windows sit ~6 months (181/184 days) after their predecessor; this
# tolerance accepts that spread while rejecting a +12-month jump when the true
# +6-month successor is absent from config["datasets"].
MIN_GAP_DAYS = 150
MAX_GAP_DAYS = 210


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="e.g. sarscov2_clades")
    parser.add_argument("--mlr-dir", default="mlr-estimates")
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


def build_pairs(windows):
    """Consecutive (W, W+1) where W+1's final date is ~6 months after W's."""
    pairs = []
    for cur, nxt in zip(windows, windows[1:]):
        if MIN_GAP_DAYS <= day_gap(nxt[2], cur[2]) <= MAX_GAP_DAYS:
            pairs.append((cur, nxt))
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
    """W+1 weekly_raw_freq at date d, folded onto W's variant set S.

    Variants absent from S are folded into "other" (when S has it); their mass
    still counts toward the denominator either way, so the residual of clades
    that emerged after W is reflected as error. Returns None when d has no
    non-null truth (smoothing-edge gaps).
    """
    s = set(variant_set)
    folded = {v: 0.0 for v in variant_set}
    total = 0.0
    seen = False
    for variant, series in truth_next.items():
        value = series.get(d)
        if value is None:
            continue
        seen = True
        total += value
        if variant in s:
            folded[variant] += value
        elif "other" in s:
            folded["other"] += value
    if not seen or total <= 0:
        return None
    return {v: folded[v] / total for v in variant_set}


def forecast_pair(cur, nxt, tau_fn, epsilon, hindcast_days, forecast_days):
    """Rows (date, lead, n, model, abs_error) for one (W, W+1) pair."""
    _, m, T = cur
    _, nm, _ = nxt
    variant_set = m["metadata"]["variants"]
    n = len(variant_set)

    ga = ff_io.variant_growth_advantages(m)          # missing / "other" -> 1.0
    modeled = ff_io.variant_modeled_frequencies(m)   # {v: {date: freq}} (median)
    truth_next = ff_io.variant_weekly_frequencies(nm)

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
    """MAE(lead) curve: mean abs_error per model within weekly lead bins."""
    bucket = defaultdict(lambda: defaultdict(list))  # bin_index -> model -> [ae]
    for _pair, _win, _nxt, _d, lead, _n, model, ae in all_rows:
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


def region_means(all_rows):
    """Per-model mean AE over the hindcast (lead <= 0) and forecast (lead > 0)."""
    acc = {"hindcast": defaultdict(list), "forecast": defaultdict(list)}
    for _pair, _win, _nxt, _d, lead, _n, model, ae in all_rows:
        acc["forecast" if lead > 0 else "hindcast"][model].append(ae)
    out = {}
    for region, by_model in acc.items():
        out[region] = {
            model: (sum(v) / len(v) if v else None)
            for model, v in {"mlr": by_model.get("mlr", []), "naive": by_model.get("naive", [])}.items()
        }
    return out


def main():
    args = parse_args()
    tau_fn = make_tau_fn(args)

    windows = load_windows(args.mlr_dir, args.dataset)
    pairs = build_pairs(windows)
    ff_io.log(
        f"{args.dataset}: {len(windows)} windows, {len(pairs)} forecastable pairs"
    )

    all_rows = []  # (pair_label, window, next_window, date, lead, n, model, ae)
    for cur, nxt in pairs:
        label = f"{cur[0]}->{nxt[0]}"
        for d, lead, n, model, ae in forecast_pair(
            cur, nxt, tau_fn, args.epsilon, args.hindcast_days, args.forecast_days
        ):
            all_rows.append((label, cur[0], nxt[0], d, lead, n, model, ae))
        ff_io.log(f"  {label}: evaluated")

    with open(args.detail_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["window", "next_window", "date", "lead_days", "n_variants", "model", "abs_error"]
        )
        for _label, window, next_window, d, lead, n, model, ae in all_rows:
            writer.writerow([window, next_window, d, lead, n, model, f"{ae:.6f}"])
    ff_io.log(f"Wrote {len(all_rows)} detail rows to {args.detail_output}")

    summary = {
        "dataset": args.dataset,
        "epsilon": args.epsilon,
        "lead_bin_days": args.lead_bin_days,
        "hindcast_days": args.hindcast_days,
        "forecast_days": args.forecast_days,
        "n_pairs": len(pairs),
        "pairs": [f"{cur[0]}->{nxt[0]}" for cur, nxt in pairs],
        "reference": 0.05,
        "overall": region_means(all_rows),
        "curve": aggregate(all_rows, args.lead_bin_days),
    }
    with open(args.summary_output, "w") as handle:
        json.dump(summary, handle, indent=2)
    ff_io.log(f"Wrote forecast-accuracy summary to {args.summary_output}")


if __name__ == "__main__":
    main()
