#!/usr/bin/env python3
"""Similarity-aware forecast scoring: error in HA1 mutation-profile space.

The label-based score (forecast_accuracy.py) treats clade labels as
exchangeable, so forecasting J.2 when J.2.4.2 sweeps is scored exactly like
forecasting B — although J.2 differs from J.2.4.2 at 7 HA1 amino acids and B at
31. For vaccine strain selection that is the wrong accounting.

This scorer keeps the predictions identical (forecast_shared.window_predictions)
and changes only the comparison. Each clade is represented by the reconstructed
HA1 amino-acid sequence at its MRCA node (clade_distances.py), a frequency
distribution over clades becomes the population's mutation profile — for each
HA1 variant position, the fraction of viruses carrying it — and the error is the
distance between the forecast profile and the observed profile.

Near-misses become cheap and real misses expensive, automatically, with no new
scoring rule: this is a fixed linear change of coordinates in front of the same
comparison.

One global map, nothing remapped
--------------------------------
The label-based score restricts truth to the *forecast* window's clade set and
renormalizes, so any clade that emerged after the date of estimation is dropped
outright — the very case this analysis is about. Here forecast and truth are each
renormalized over their own window's named clades and then both are embedded in
one global profile space spanning every clade in the tree. Their supports never
have to agree. Nothing is remapped, dropped, or relabeled.

Units
-----
Error is reported in amino acids: the distance between two pure clades is the
count of HA1 positions at which their MRCA sequences differ. Reported results use
the raw HA1 distance (`--eps 0`), alongside the label-based comparator `mae`,
the mean absolute frequency error. Because the forecast and the truth sit on
different clade sets, that mean is taken over the UNION of the two sets; the
divisor is the same for both models, so it cannot affect the comparison between
them, only the level. `tv` (total variation, half the summed absolute error) is
carried alongside because it is what the verification limit below converges to.

Note that the raw distance deliberately does not separate clades whose MRCA
sequences are identical (H3N2 has one such pair, J and J.3): to this score they
are the same virus, which is the point of scoring HA1 rather than labels.

`eps` exists for verification, not for reporting. It adds a floor of separation,
in amino acids, between every pair of distinct clades via one identity column per
clade. As eps grows, clade identity is all that survives and the normalized error
converges to `tv` under the default L1 norm (to the Euclidean difference norm
under L2, which is why L1 is the default). `--self-test` asserts exactly that, so
a mistake in the encoding, the scaling or the normalization cannot pass silently.
"""

import argparse
import csv
import json
from collections import defaultdict

import numpy as np

import ff_io
from forecast_shared import (
    build_pairs,
    load_windows,
    make_tau_fn,
    named_variants,
    window_predictions,
)

# Clades with no reconstructable HA1 sequence, so no position in profile space.
# "other" is the collapsed rare-clade bucket, already excluded by the label-based
# score. "unassigned" is the tree's catch-all for tips outside any named clade;
# the label-based score does treat it as a scoreable clade, but it has no MRCA
# and therefore no profile, so this score drops it and renormalizes without it.
# In H3N2 it averages 0.35% of sequences (above 1% in 5.5% of region-weeks).
UNSCOREABLE = {"other", "unassigned"}


def scoreable(freqs):
    """Drop clades with no HA1 profile and renormalize the rest to sum 1."""
    kept = {c: v for c, v in freqs.items() if c not in UNSCOREABLE}
    total = sum(kept.values())
    if not kept or total <= 0:
        return None
    return {c: v / total for c, v in kept.items()}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--dataset", required=True, help="e.g. h3n2_clades")
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument("--profile-table", required=True, help="clade x column 0/1 table")
    parser.add_argument("--distance-matrix", required=True, help="pairwise clade aa distances")
    parser.add_argument(
        "--eps",
        type=float,
        action="append",
        default=None,
        help="verification knob: amino-acid floor between distinct clades. "
        "Reported results use 0; repeatable (default: 0)",
    )
    parser.add_argument("--norm", choices=["l1", "l2"], default="l1")
    parser.add_argument("--slide-days", type=int, default=180)
    parser.add_argument("--generation-time", type=float, default=3.2)
    parser.add_argument("--generation-time-pre-omicron", type=float, default=None)
    parser.add_argument("--variant-classification", choices=["clades", "lineages"], default=None)
    parser.add_argument("--aliasing", type=str, default=None)
    parser.add_argument("--epsilon", type=float, default=1e-3)
    parser.add_argument(
        "--endpoint",
        type=int,
        action="append",
        default=None,
        help="lead times to summarize as horizon endpoints (default: 90, 180, 365)",
    )
    parser.add_argument(
        "--endpoint-span",
        type=int,
        default=6,
        help="width in days of each endpoint band, ending at the endpoint",
    )
    parser.add_argument("--lead-bin-days", type=int, default=7)
    parser.add_argument("--hindcast-days", type=int, default=90)
    parser.add_argument("--forecast-days", type=int, default=180)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run the scoring-transform correctness checks, then exit",
    )
    parser.add_argument("--detail-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def load_profile_table(path):
    """(clade index, column names, 0/1 matrix) from clade_distances.py."""
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        columns = next(reader)[1:]
        clades, rows = [], []
        for record in reader:
            clades.append(record[0])
            rows.append([int(value) for value in record[1:]])
    return {clade: i for i, clade in enumerate(clades)}, columns, np.array(rows, dtype=float)


def load_distances(path):
    """{clade_a: {clade_b: aa distance}} from clade_distances.py."""
    table = defaultdict(dict)
    with open(path) as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            table[record["clade_a"]][record["clade_b"]] = int(record["aa_dist"])
    return table


def load_max_distance(path):
    """Largest pairwise clade distance, which sets the normalization scale."""
    worst = 0
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for record in reader:
            worst = max(worst, int(record["aa_dist"]))
    return worst


def augmented_matrix(profile, eps, norm):
    """Profile columns scaled to amino acids, plus one identity column per clade.

    Two clades differ in exactly two profile columns per differing position, so
    the scale is chosen to make the distance between two pure clades come out as
    the plain amino-acid count under the chosen norm:

      L1: columns * 1/2,        identity eps/2        -> distance = d + eps
      L2: columns * 1/sqrt(2),  identity sqrt(eps/2)  -> distance = sqrt(d + eps)

    The identity block is what keeps the map injective, so distinct clade
    mixtures can never collapse onto the same profile.
    """
    n_clades = profile.shape[0]
    if norm == "l1":
        scaled = profile * 0.5
        identity = np.eye(n_clades) * (eps / 2.0)
    else:
        scaled = profile / np.sqrt(2.0)
        identity = np.eye(n_clades) * np.sqrt(eps / 2.0)
    return np.hstack([scaled, identity])


def worst_case(max_distance, eps, norm):
    """Distance between the two most distant clades, for normalizing to [0, 1]."""
    return (max_distance + eps) if norm == "l1" else np.sqrt(max_distance + eps)


def to_vector(freqs, index, size):
    """{clade: frequency} -> dense vector over the global clade index.

    Mass on a clade absent from the map (there should be none once the distance
    matrix is built from the same MLR windows) is reported by the caller rather
    than silently dropped.
    """
    vector = np.zeros(size)
    missing = 0.0
    for clade, value in freqs.items():
        position = index.get(clade)
        if position is None:
            missing += value
        else:
            vector[position] += value
    return vector, missing


def normalize_over(series, variants, d):
    """W's or W+1's named-clade frequencies at date d, renormalized to sum 1.

    Each side is renormalized over *its own* window's named clades. The two
    supports need not agree — that is the point of the global map.
    """
    kept = {}
    for variant in variants:
        value = series.get(variant, {}).get(d)
        if value is not None:
            kept[variant] = value
    total = sum(kept.values())
    if not kept or total <= 0:
        return None
    return {variant: value / total for variant, value in kept.items()}


def score(delta, matrices, worst, norm):
    """{eps: (error in amino acids, error normalized to [0, 1])} for one difference vector."""
    out = {}
    for eps, matrix in matrices.items():
        projected = delta @ matrix
        error = np.abs(projected).sum() if norm == "l1" else np.sqrt((projected ** 2).sum())
        out[eps] = (float(error), float(error / worst[eps]))
    return out


def aggregate(rows, bin_days, eps_values):
    """Mean error per model within lead bins, pooled over regions and window pairs."""
    bucket = defaultdict(lambda: defaultdict(list))
    for (_region, _pair, _win, _nxt, _d, lead, _nf, _nt, _nu,
         eps, model, aa, normed, mae, tv) in rows:
        bucket[lead // bin_days][(model, eps)].append(aa)
        bucket[lead // bin_days][(model, "norm", eps)].append(normed)
        bucket[lead // bin_days][(model, "mae")].append(mae)
        bucket[lead // bin_days][(model, "tv")].append(tv)
    curve = []
    for b in sorted(bucket):
        entry = {
            "lead": b * bin_days + (bin_days - 1) / 2.0,
            "lead_lo": b * bin_days,
            "lead_hi": b * bin_days + bin_days - 1,
        }
        for model in ("mlr", "naive"):
            for eps in eps_values:
                for key, values in (
                    (f"{model}_eps{eps:g}", bucket[b].get((model, eps), [])),
                    (f"{model}_norm_eps{eps:g}", bucket[b].get((model, "norm", eps), [])),
                ):
                    entry[key] = sum(values) / len(values) if values else None
            for key in ("mae", "tv"):
                values = bucket[b].get((model, key), [])
                entry[f"{model}_{key}"] = sum(values) / len(values) if values else None
            entry[f"n_{model}"] = len(bucket[b].get((model, eps_values[0]), []))
        curve.append(entry)
    return curve


def horizon_means(rows, eps_values):
    """Per-model mean error over the hindcast (lead <= 0) and forecast (lead > 0)."""
    acc = {"hindcast": defaultdict(list), "forecast": defaultdict(list)}
    for (_region, _pair, _win, _nxt, _d, lead, _nf, _nt, _nu,
         eps, model, aa, normed, mae, tv) in rows:
        horizon = "forecast" if lead > 0 else "hindcast"
        acc[horizon][(model, eps)].append(aa)
        acc[horizon][(model, "norm", eps)].append(normed)
        acc[horizon][(model, "mae")].append(mae)
        acc[horizon][(model, "tv")].append(tv)
    out = {}
    for horizon, by_key in acc.items():
        entry = {}
        for model in ("mlr", "naive"):
            for eps in eps_values:
                for key, values in (
                    (f"{model}_eps{eps:g}", by_key.get((model, eps), [])),
                    (f"{model}_norm_eps{eps:g}", by_key.get((model, "norm", eps), [])),
                ):
                    entry[key] = sum(values) / len(values) if values else None
            for key in ("mae", "tv"):
                values = by_key.get((model, key), [])
                entry[f"{model}_{key}"] = sum(values) / len(values) if values else None
        out[horizon] = entry
    return out


def label_limit(delta, norm):
    """The label-based score each norm collapses to once clades are equidistant.

    As eps grows the identity block swamps the profile columns and clade identity
    is all that survives. What remains differs by norm, and the difference is
    real rather than an artifact:

      L1 -> total variation, 0.5*sum|delta|, which is the label-based comparator
            reported as `tv` (and, up to the 1/n factor, the published MAE)
      L2 -> the Euclidean norm of the difference, ||delta||/sqrt(2)

    L1 is the default precisely because its limit is the metric the paper
    already reports.
    """
    if norm == "l1":
        return 0.5 * float(np.abs(delta).sum())
    return float(np.sqrt((delta ** 2).sum()) / np.sqrt(2.0))


def ordering_test(profile, index, distances, norm):
    """Shifting forecast mass to a more distant clade must cost more.

    The brief's sanity check, kept as a permanent test: hold the winner fixed and
    move the forecast onto its nearest relative, then onto a mid-distance clade,
    then onto the most distant one. The scores must come out in that order.
    """
    size = profile.shape[0]
    matrix = augmented_matrix(profile, 0.0, norm)

    winner = max(distances, key=lambda c: len(distances[c]))
    ranked = sorted((c for c in distances[winner] if c != winner), key=lambda c: distances[winner][c])
    picks = [ranked[0], ranked[len(ranked) // 2], ranked[-1]]

    truth = np.zeros(size)
    truth[index[winner]] = 1.0
    scores = []
    for clade in picks:
        forecast = np.zeros(size)
        forecast[index[clade]] = 0.8
        forecast[index[winner]] = 0.2
        projected = (forecast - truth) @ matrix
        error = np.abs(projected).sum() if norm == "l1" else np.sqrt((projected ** 2).sum())
        scores.append((clade, distances[winner][clade], float(error)))

    if not all(scores[i][2] < scores[i + 1][2] for i in range(len(scores) - 1)):
        raise SystemExit(f"SELF-TEST FAILED ({norm}): scores not ordered by distance: {scores}")
    detail = ", ".join(f"{c} ({d} aa) -> {e:.2f}" for c, d, e in scores)
    ff_io.log(f"self-test passed ({norm}): truth {winner}, 80% of mass on {detail}")


def endpoint_stats(rows, targets, span, forecast_days):
    """Per-horizon mean error and win rate, from unrounded in-memory values.

    The win rate is decided by strict inequality, so it must be computed here
    rather than re-derived from the written TSV: at long lead the two models'
    errors agree to several decimals often enough that the file's rounding moves
    the tie count by whole percentage points.

    Each endpoint covers leads in (target - span, target].
    """
    out = []
    eps0 = min(r[9] for r in rows) if rows else None
    for target in targets:
        if target > forecast_days:
            continue
        by = defaultdict(dict)
        for (region, _pair, win, nxt, d, lead, _nf, _nt, _nu,
             eps, model, aa, _normed, mae, _tv) in rows:
            if eps != eps0 or not (target - span < lead <= target):
                continue
            by[(region, win, nxt, d)][model] = (aa, mae)
        paired = [v for v in by.values() if "mlr" in v and "naive" in v]
        if not paired:
            continue
        entry = {"lead": target, "lead_lo": target - span + 1, "n": len(paired)}
        for i, key in ((0, "aa"), (1, "mae")):
            mlr = sum(v["mlr"][i] for v in paired) / len(paired)
            naive = sum(v["naive"][i] for v in paired) / len(paired)
            wins = sum(1 for v in paired if v["mlr"][i] < v["naive"][i])
            ties = sum(1 for v in paired if v["mlr"][i] == v["naive"][i])
            entry[key] = {
                "mlr": mlr,
                "naive": naive,
                "advantage_pct": (naive - mlr) / naive * 100 if naive else None,
                "win_rate_pct": wins / len(paired) * 100,
                "tie_rate_pct": ties / len(paired) * 100,
            }
        out.append(entry)
    return out


def self_test(profile, max_distance, norm):
    """The large-eps limit must reproduce this norm's label-based score.

    This is the regression test for the whole transform: if the identity block,
    the scaling, or the normalization is wrong, the limit will not land.
    """
    rng = np.random.default_rng(0)
    n = profile.shape[0]
    failures = 0
    for _ in range(200):
        p = rng.dirichlet(np.ones(n) * 0.3)
        q = rng.dirichlet(np.ones(n) * 0.3)
        target = label_limit(p - q, norm)
        eps = 1e6
        matrix = augmented_matrix(profile, eps, norm)
        worst = worst_case(max_distance, eps, norm)
        projected = (p - q) @ matrix
        error = np.abs(projected).sum() if norm == "l1" else np.sqrt((projected ** 2).sum())
        got = error / worst
        if target > 0 and abs(got - target) / target > 1e-3:
            failures += 1
    if failures:
        raise SystemExit(
            f"SELF-TEST FAILED ({norm}): large-eps limit missed the label-based "
            f"score in {failures}/200 draws"
        )
    label = "total variation" if norm == "l1" else "the Euclidean difference norm"
    ff_io.log(f"self-test passed ({norm}): large-eps limit reproduces {label} in 200/200 draws")

    # And a pure-clade pair must sit at exactly the amino-acid distance apart.
    matrix = augmented_matrix(profile, 0.0, norm)
    a, b = 0, 1
    delta = np.zeros(n)
    delta[a], delta[b] = 1.0, -1.0
    projected = delta @ matrix
    error = np.abs(projected).sum() if norm == "l1" else np.sqrt((projected ** 2).sum())
    expected = np.abs(profile[a] - profile[b]).sum() / 2.0
    if norm == "l2":
        expected = np.sqrt(expected)
    if abs(error - expected) > 1e-9:
        raise SystemExit(f"SELF-TEST FAILED: pure-clade distance {error} != {expected}")
    ff_io.log(f"self-test passed ({norm}): pure-clade distance reads in amino acids")


def main():
    args = parse_args()
    eps_values = args.eps if args.eps is not None else [0.0]
    args.endpoint = args.endpoint if args.endpoint is not None else [90, 180, 365]
    eps_values = sorted(set(eps_values))

    index, columns, profile = load_profile_table(args.profile_table)
    max_distance = load_max_distance(args.distance_matrix)
    ff_io.log(
        f"profile map: {profile.shape[0]} clades x {len(columns)} columns, "
        f"max pairwise distance {max_distance} aa"
    )

    if args.self_test:
        self_test(profile, max_distance, args.norm)
        ordering_test(profile, index, load_distances(args.distance_matrix), args.norm)
        return

    matrices = {eps: augmented_matrix(profile, eps, args.norm) for eps in eps_values}
    worst = {eps: worst_case(max_distance, eps, args.norm) for eps in eps_values}
    size = profile.shape[0]

    tau_fn = make_tau_fn(args)
    windows = load_windows(args.mlr_dir, args.dataset)
    all_regions = sorted(set().union(*[set(ff_io.regions(w[1])) for w in windows])) if windows else []

    rows = []
    pairs_by_region = {}
    dropped_mass = 0.0
    dropped_count = 0
    for region in all_regions:
        region_windows = [w for w in windows if region in ff_io.regions(w[1])]
        pairs = build_pairs(region_windows, args.slide_days)
        pairs_by_region[region] = [f"{cur[0]}->{nxt[0]}" for cur, nxt in pairs]
        for cur, nxt in pairs:
            label = f"{cur[0]}->{nxt[0]}"
            _, nm, _ = nxt
            truth_variants = [v for v in named_variants(nm) if v not in UNSCOREABLE]
            truth_series = ff_io.variant_weekly_frequencies(nm, location=region)

            for d, lead, variant_set, pred_mlr, pred_naive in window_predictions(
                cur, nxt, tau_fn, args.epsilon, args.hindcast_days, args.forecast_days,
                location=region,
            ):
                truth = normalize_over(truth_series, truth_variants, d)
                if truth is None:
                    continue
                prediction_pair = [
                    (name, scoreable(values))
                    for name, values in (("mlr", pred_mlr), ("naive", pred_naive))
                ]
                if any(values is None for _, values in prediction_pair):
                    continue
                truth_vector, missing = to_vector(truth, index, size)
                if missing:
                    dropped_mass += missing
                    dropped_count += 1
                # Mean absolute error is taken over the union of the forecast's and
                # the truth's clade sets, since the two need not agree. Both models
                # share the same support, so the divisor is common to them.
                n_union = len(set(truth) | set(prediction_pair[0][1]))
                for model, prediction in prediction_pair:
                    pred_vector, pred_missing = to_vector(prediction, index, size)
                    if pred_missing:
                        dropped_mass += pred_missing
                        dropped_count += 1
                    delta = pred_vector - truth_vector
                    total = float(np.abs(delta).sum())
                    tv = 0.5 * total
                    mae = total / n_union if n_union else None
                    for eps, (aa, normed) in score(delta, matrices, worst, args.norm).items():
                        rows.append((
                            region, label, cur[0], nxt[0], d, lead,
                            len(prediction), len(truth_variants), n_union,
                            eps, model, aa, normed, mae, tv,
                        ))
        ff_io.log(f"  {region}: {len(region_windows)} windows, {len(pairs)} pairs")

    total_pairs = sum(len(p) for p in pairs_by_region.values())
    ff_io.log(
        f"{args.dataset}: {len(all_regions)} regions, {total_pairs} region-window pairs "
        f"(slide {args.slide_days}d, horizon +{args.forecast_days}d / -{args.hindcast_days}d, "
        f"norm {args.norm}, eps {eps_values})"
    )
    if dropped_count:
        ff_io.log(
            f"WARNING: {dropped_count} vectors carried mass on clades absent from the "
            f"profile map (total {dropped_mass:.3f})"
        )

    with open(args.detail_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "region", "window", "next_window", "date", "lead_days",
            "n_forecast_clades", "n_truth_clades", "n_union_clades", "eps", "model",
            "error_aa", "error_norm", "mae", "tv",
        ])
        for (region, _label, window, next_window, d, lead, nf, nt, nu,
             eps, model, aa, normed, mae, tv) in rows:
            writer.writerow([
                region, window, next_window, d, lead, nf, nt, nu,
                # These columns span two orders of magnitude, and win rates are
                # decided by strict inequality, so they are all written at a
                # precision where near-ties cannot round together into a tie.
                f"{eps:g}", model, f"{aa:.8f}", f"{normed:.8f}", f"{mae:.10f}", f"{tv:.10f}",
            ])
    ff_io.log(f"Wrote {len(rows)} detail rows to {args.detail_output}")

    summary = {
        "dataset": args.dataset,
        "norm": args.norm,
        "eps_values": eps_values,
        "max_distance": max_distance,
        "n_clades": size,
        "n_profile_columns": len(columns),
        "slide_days": args.slide_days,
        "lead_bin_days": args.lead_bin_days,
        "hindcast_days": args.hindcast_days,
        "forecast_days": args.forecast_days,
        "regions": all_regions,
        "n_pairs": total_pairs,
        "pairs_by_region": pairs_by_region,
        "overall": horizon_means(rows, eps_values),
        "endpoints": endpoint_stats(rows, args.endpoint, args.endpoint_span, args.forecast_days),
        "curve": aggregate(rows, args.lead_bin_days, eps_values),
    }
    with open(args.summary_output, "w") as handle:
        json.dump(summary, handle, indent=2)
    ff_io.log(f"Wrote similarity summary to {args.summary_output}")


if __name__ == "__main__":
    main()
