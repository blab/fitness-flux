#!/usr/bin/env python3
"""Window pairing and forward projection shared by the forecast scorers.

Both forecast_accuracy.py (label-based mean absolute error) and
forecast_similarity.py (HA1 similarity-aware error) make the *same* predictions
and differ only in how they compare those predictions to truth. That common part
lives here so the two scores can never drift apart.

Predictions, over a window's *named* clades (the "other" bucket excluded and the
rest renormalized to sum 1):

  MLR:   lead <= 0  -> the in-window modeled freq_v(t) (site "freq", median)
         lead >  0  -> softmax_v[ ln max(freq_v(T), eps) + delta_v*lead ]
  Naive: lead <= 0  -> identical to MLR (both are the in-window fit)
         lead >  0  -> freq_v(T) held constant (flat forward from T)

delta_v (per-day log-growth slope) is recovered exactly from the exported growth
advantage: the MLR export writes ga_v = exp(delta_v * tau_v) with the same
per-variant tau_v used at fit time (scripts/generation_time.py), so
delta_v = ln(ga_v) / tau_v. Since evofr's freq_v(t) = softmax_v(intercept +
slope*t) and softmax is shift-invariant, softmax_v[ln freq_v(T) + delta_v*lead]
reconstructs the model's own forward projection.
"""

import math
import os
import sys
from datetime import date

import ff_io

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from generation_time import generation_time_for, load_aliasor  # noqa: E402

# A window's successor sits ~one slide later. We accept a successor whose final
# data date is within +/-20% of the slide, which spans the real spread (e.g.
# 181-184 d for a 6-month slide, 324-407 d for a 1-year slide) while rejecting a
# doubled jump when the true successor is absent, and same-year duplicates from
# bridge windows.
GAP_TOLERANCE = 0.20


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


def named_variants(mlr):
    """A window's named clades: metadata variants minus the "other" bucket."""
    return [v for v in mlr["metadata"]["variants"] if v != "other"]


def window_predictions(cur, nxt, tau_fn, epsilon, hindcast_days, forecast_days, location=None):
    """Yield (date, lead, variant_set, pred_mlr, pred_naive) for one (W, W+1) pair.

    Predictions only — truth and error are the caller's business, which is what
    lets the label-based and similarity-aware scorers share this code. Yields
    nothing when W has fewer than two named clades (relative frequencies are
    undefined). Dates come from W+1's grid, filtered to the lead window.
    """
    _, m, T = cur
    _, nm, _ = nxt
    variant_set = named_variants(m)
    if len(variant_set) < 2:
        return

    ga = ff_io.variant_growth_advantages(m, location=location)          # missing / "other" -> 1.0
    modeled = ff_io.variant_modeled_frequencies(m, location=location)   # {v: {date: freq}} (median)

    freqT = normalize({v: (modeled.get(v, {}).get(T) or 0.0) for v in variant_set})
    delta = {v: math.log(max(ga.get(v, 1.0), epsilon)) / tau_fn(v) for v in variant_set}
    logf0 = {v: math.log(max(freqT[v], epsilon)) for v in variant_set}

    for d in nm["metadata"]["dates"]:
        lead = day_gap(d, T)
        if lead < -hindcast_days or lead > forecast_days:
            continue
        if lead <= 0:
            # Hindcast: both models are the in-window fit, so they coincide.
            pred_mlr = normalize({v: (modeled.get(v, {}).get(d) or 0.0) for v in variant_set})
            pred_naive = pred_mlr
        else:
            pred_mlr = softmax({v: logf0[v] + delta[v] * lead for v in variant_set})
            pred_naive = freqT
        yield d, lead, variant_set, pred_mlr, pred_naive
