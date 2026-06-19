#!/usr/bin/env python3
"""Fit linear models relating branch mutation counts to fitness change.

Port of the "Linear model across predictors" and "Slope through time" sections of
``lineage-deltas.nb``. Two outputs:

  * linear_model_coefficients.tsv -- multiple regression of delta log fitness on
    the per-region mutation-count deltas (ORF1ab, NTD, RBD, S2, accessory), fit
    over all branches and separately over early vs late seasons.
  * slope_through_time.tsv -- per-season simple-regression slope of delta log
    fitness on spike and non-spike mutation deltas.
"""

import argparse
import csv

import numpy as np

import ld_io

# regression design: delta_log_fitness ~ these mutation-region deltas
PREDICTOR_COLUMNS = [
    "delta_orf1ab_muts",
    "delta_ntd_muts",
    "delta_rbd_muts",
    "delta_s2_muts",
    "delta_accessory_muts",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--branch-deltas", required=True)
    parser.add_argument("--coefficients-output", required=True)
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


def fit_multiple(rows):
    """Return (coefficients dict incl. intercept, r_squared, n)."""
    if len(rows) <= len(PREDICTOR_COLUMNS) + 1:
        return None
    x = np.array(
        [[1.0] + [float(r[c]) for c in PREDICTOR_COLUMNS] for r in rows]
    )
    y = np.array([float(r["delta_log_fitness"]) for r in rows])
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    residuals = y - x @ coef
    ss_res = float(np.sum(residuals**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    terms = ["intercept"] + PREDICTOR_COLUMNS
    return dict(zip(terms, coef)), r_squared, len(rows)


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

    with open(args.coefficients_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["group", "term", "coefficient", "r_squared", "n"])
        for group, group_rows in groups.items():
            fit = fit_multiple(group_rows)
            if fit is None:
                continue
            coefficients, r_squared, n = fit
            for term, value in coefficients.items():
                writer.writerow(
                    [group, term, f"{value:.6f}", f"{r_squared:.4f}", n]
                )
    ld_io.log(f"Wrote linear-model coefficients to {args.coefficients_output}")

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
