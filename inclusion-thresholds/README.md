# Inclusion thresholds for MLR variant modeling

How to set the per-window thresholds that decide which variants are modeled
separately vs collapsed before MLR fitting. Two write-ups, with figures:

- [`clade-analysis.md`](clade-analysis.md) — Nextstrain **clades** (`clade_min_count` /
  `clade_min_freq`), which collapse rare clades into a universal `other`.
- [`lineage-analysis.md`](lineage-analysis.md) — Pango **lineages** (`collapse_threshold` +
  `min_core_count`), which roll rare lineages up into their parent.

The two reach **opposite** conclusions, for a principled reason. In both, scoring
the fit on its own labels (cross-entropy, or equally MSE/MAE) is confounded — it
rewards collapsing — so neither selects on it.

- **Clades**: the binding failure mode is over-collapse into an incoherent `other`
  bucket, which scales with **frequency**. → use a mean-frequency threshold,
  `clade_min_freq = 0.1 %` (implemented as `--clade-min-freq` in `prepare-data.py`,
  set for all clade datasets in `defaults/config.yaml`).
- **Lineages**: parent-rollup makes over-collapse gentle, so the binding mode is
  estimation **noise**, which scales with absolute **count** (growth-rate
  SE ≈ 0.084/√count). → keep `collapse_threshold` a count (`500`; ~200–500 supported,
  a frequency would be worse). A second, structural failure mode — *relic residuals*
  (near-extinct ancestors padded out by scattered descendants) that corrupt the
  lineage mutation→fitness deltas — is fixed by a self-supporting `min_core_count`
  floor (`200`, `--min-core-count` in `collapse-lineage-counts.py`). Both are tuned
  for the lineage **deltas** analysis; clades carry the flux.

## Layout

```
scripts/sweep_clade_min_seq.py        # for each (dataset, theta): prepare-data (+ collapse for
                                      #   lineages) + MLR -> scratch/; clades use MAP, lineages MAP
scripts/mlr-config-map.yaml           # MAP config for the clade sweep
scripts/mlr-config-lineage-map.yaml   # MAP config for the lineage sweep (pivot auto)
scripts/mlr-config-lineage-nuts.yaml  # NUTS config for lineage re-validation at the recommended theta
scripts/score_sweep.py                # clades: coarse/fixed CE + per-clade peaks -> results/
scripts/score_lineage_sweep.py        # lineages: coarse/fixed CE, analytic SE noise, fitness-signal-lost
scripts/plot_sweep.py                 # clade figures + recommendation
scripts/plot_lineage_sweep.py         # lineage figures + recommendation
scripts/verify_offramp.py             # tunes min_core_count to the lineage deltas correlation (12 seasons)
results/                              # *.tsv metrics (not versioned)
figures/                              # ce_confound/overcollapse/clade_peaks + lineage_* (versioned)
scratch/                              # per-threshold counts + MLR fits + cached alias_key.json (not versioned)
```

## Reproduce

From the repository root (reuses the cached, threshold-independent
`sequence-counts/{dataset}/seq_counts.tsv`; does not touch the production config or
canonical outputs):

```
# clades
python inclusion-thresholds/scripts/sweep_clade_min_seq.py
python inclusion-thresholds/scripts/score_sweep.py
python inclusion-thresholds/scripts/plot_sweep.py

# lineages (per window: sarscov2_lineages_{2020-21,2021-22,2023,2025})
python inclusion-thresholds/scripts/sweep_clade_min_seq.py --dataset sarscov2_lineages_2021-22
python inclusion-thresholds/scripts/score_lineage_sweep.py
python inclusion-thresholds/scripts/plot_lineage_sweep.py
```

The sweep is idempotent — it skips any `(dataset, theta)` whose fit already exists.
Datasets and the count/frequency grids are defined at the top of
`sweep_clade_min_seq.py`. The lineage collapse step reads a cached
`scratch/alias_key.json` (pre-fetched once) via `collapse-lineage-counts.py --aliasing`.
