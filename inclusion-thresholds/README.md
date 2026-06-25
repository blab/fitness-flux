# Principled selection of the rare-clade collapsing threshold

A standalone study of how to set the threshold below which a clade is collapsed
into `other` before MLR fitting — the per-window count lever `clade_min_count`
(renamed from `clade_min_seq`) and its frequency-based replacement
`clade_min_freq`. The write-up, with figures, is in [`clade-analysis.md`](clade-analysis.md).

Short version: scoring the fit (cross-entropy, or equally MSE/MAE) on each fit's
own labels is confounded — it rewards collapsing. The robust criterion is the
**peak smoothed frequency of the largest clade folded into `other`** (boundary
between real clades and artifacts ≈ 1 %). Because a clade's peak tracks its mean
frequency tightly (r = 0.94), the simple raw-count lever is a **mean-frequency
threshold**: this analysis selects `clade_min_freq = 0.1 %`, now implemented as
`--clade-min-freq` in `prepare-data.py` and set for all clade datasets in
`defaults/config.yaml`.

## Layout

```
scripts/sweep_clade_min_seq.py   # for each (dataset, theta): prepare-data + MAP MLR -> scratch/
scripts/mlr-config-map.yaml      # fast MAP config for the sweep (production uses NUTS)
scripts/score_sweep.py           # coarse & fixed-partition CE/BIC + per-clade peaks -> results/
scripts/plot_sweep.py            # figures/ + recommendation table
results/sweep_metrics.tsv        # per (dataset, theta) metrics
results/clade_peaks.tsv          # per-clade window-total, mean freq, smoothed peak
figures/                         # ce_confound, overcollapse, clade_peaks
scratch/                         # per-threshold prepared counts + MLR fits (not versioned)
```

## Reproduce

From the repository root (reuses the cached, threshold-independent
`sequence-counts/{dataset}/seq_counts.tsv`; does not touch the production config
or canonical outputs):

```
python inclusion-thresholds/scripts/sweep_clade_min_seq.py
python inclusion-thresholds/scripts/score_sweep.py
python inclusion-thresholds/scripts/plot_sweep.py
```

The sweep is idempotent — it skips any `(dataset, theta)` whose fit already
exists, so it can be resumed or extended. Datasets and the count/frequency grids
are defined at the top of `sweep_clade_min_seq.py`.
