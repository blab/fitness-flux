#!/usr/bin/env python3
"""Fit linear models relating branch mutation counts to fitness change.

Port of the "Linear model across predictors" and "Slope through time" sections of
``lineage-deltas.nb``. Three outputs:

  * linear_model_coefficients.tsv -- multiple regression of delta log fitness on
    the change in substitution count across four non-overlapping genome regions
    (spike RBD, spike S1 outside the RBD, ORF1ab, accessory), with per-term
    estimate / standard error / t / p, fit over all branches and separately over
    early vs late seasons. The multiple regression isolates each region's partial
    contribution, removing the confounding present in the marginal scatters (a
    more-evolved lineage accrues substitutions everywhere at once).
  * linear_model_predictions.tsv -- per-branch fitted value from the pooled
    ("all") model, for the predicted-vs-observed scatter.
  * slope_through_time.tsv -- per-season simple-regression slope of delta log
    fitness on spike and non-spike mutation deltas.
"""

import argparse
import csv

import numpy as np
from scipy import stats

import ld_io

# regression design: delta_log_fitness ~ these non-overlapping region deltas
PREDICTOR_COLUMNS = [
    "delta_rbd_muts",
    "delta_s1_nonrbd_muts",
    "delta_orf1ab_muts",
    "delta_accessory_muts",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--branch-deltas", required=True)
    parser.add_argument("--coefficients-output", required=True)
    parser.add_argument("--predictions-output", required=True)
    parser.add_argument("--slope-output", required=True)
    parser.add_argument(
        "--late-from",
        default="2022",
        help="seasons whose label sorts >= this are 'late' (else 'early')",
    )
    return parser.parse_args()


def read_rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def design_matrix(rows):
    return np.array(
        [[1.0] + [float(r[c]) for c in PREDICTOR_COLUMNS] for r in rows]
    )


def fit_multiple(rows):
    """Fit OLS ``delta_log_fitness ~ 1 + PREDICTOR_COLUMNS``.

    Return a dict with per-term estimate / standard error / t / p (ordinary OLS
    inference: coefficient covariance ``sigma^2 (X'X)^-1``, two-sided Student-t
    p-values on ``n - k`` degrees of freedom), plus ``r_squared``, ``n`` and the
    raw coefficient vector for generating predictions. ``None`` if too few rows.
    """
    k = len(PREDICTOR_COLUMNS) + 1  # predictors + intercept
    n = len(rows)
    if n <= k:
        return None
    x = design_matrix(rows)
    y = np.array([float(r["delta_log_fitness"]) for r in rows])
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    residuals = y - x @ coef
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    dof = n - k
    sigma2 = ss_res / dof
    cov = sigma2 * np.linalg.inv(x.T @ x)
    se = np.sqrt(np.diag(cov))
    tvals = coef / se
    pvals = 2.0 * stats.t.sf(np.abs(tvals), dof)

    terms = ["intercept"] + PREDICTOR_COLUMNS
    term_rows = [
        {
            "term": term,
            "estimate": float(coef[i]),
            "se": float(se[i]),
            "t": float(tvals[i]),
            "p": float(pvals[i]),
        }
        for i, term in enumerate(terms)
    ]
    return {"terms": term_rows, "r_squared": r_squared, "n": n, "coef": coef}


def simple_slope(rows, column):
    if len(rows) < 3:
        return None
    x = np.array([float(r[column]) for r in rows])
    y = np.array([float(r["delta_log_fitness"]) for r in rows])
    if np.ptp(x) == 0:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope)


def main():
    args = parse_args()
    rows = read_rows(args.branch_deltas)

    groups = {
        "all": rows,
        "early": [r for r in rows if r["timepoint"] < args.late_from],
        "late": [r for r in rows if r["timepoint"] >= args.late_from],
    }

    fits = {group: fit_multiple(group_rows) for group, group_rows in groups.items()}

    with open(args.coefficients_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["group", "term", "estimate", "se", "t", "p", "r_squared", "n"]
        )
        for group, fit in fits.items():
            if fit is None:
                continue
            for term in fit["terms"]:
                writer.writerow(
                    [
                        group,
                        term["term"],
                        f"{term['estimate']:.6f}",
                        f"{term['se']:.6f}",
                        f"{term['t']:.4f}",
                        f"{term['p']:.3e}",
                        f"{fit['r_squared']:.4f}",
                        fit["n"],
                    ]
                )
    ld_io.log(f"Wrote linear-model coefficients to {args.coefficients_output}")

    # Per-branch fitted values from the pooled ("all") model.
    all_fit = fits["all"]
    with open(args.predictions_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["timepoint", "parent", "child", "lm_predicted"])
        if all_fit is not None:
            predicted = design_matrix(rows) @ all_fit["coef"]
            for row, pred in zip(rows, predicted):
                writer.writerow(
                    [row["timepoint"], row["parent"], row["child"], f"{pred:.6f}"]
                )
    ld_io.log(f"Wrote linear-model predictions to {args.predictions_output}")

    timepoints = sorted({r["timepoint"] for r in rows})
    with open(args.slope_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["timepoint", "predictor", "slope", "n"])
        for timepoint in timepoints:
            season_rows = [r for r in rows if r["timepoint"] == timepoint]
            for column, label in [
                ("delta_spike_muts", "spike"),
                ("delta_nonspike_muts", "nonspike"),
                ("delta_nuc_muts", "all"),
            ]:
                slope = simple_slope(season_rows, column)
                if slope is not None:
                    writer.writerow(
                        [timepoint, label, f"{slope:.6f}", len(season_rows)]
                    )
    ld_io.log(f"Wrote slope-through-time to {args.slope_output}")


if __name__ == "__main__":
    main()
