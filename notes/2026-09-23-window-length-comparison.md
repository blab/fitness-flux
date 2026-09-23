# H3N2 window length: 1-year vs 2-year, scored by HA1 distance

**Verdict: keep 2-year windows.** No difference in forecast skill was detectable
between 1-year and 2-year H3N2 MLR windows. Point estimates favour 2-year at 90
and 180 days, but every comparison's confidence interval straddles zero once
uncertainty is clustered properly. 2-year stays on the status quo, and on having
roughly twice the sequences per fit.

The production config was never changed; the 1-year arm ran from a config overlay
and has been removed. This note exists so the next person does not have to redo
the work — the previous 1-year experiment left no config, branch or note, only a
sentence in the commit message of `4f8461f`.

## What was run

Production 2-year windows against a standalone 1-year configuration under the
same thresholds (`clade_min_count: 50`, `clade_min_freq: 0.001`,
`location_min_seq: 500`), hierarchical multi-region MLR with `pool_scale: 0.25`,
33 windows spanning 2017-01 to 2025-10. Scored with
`fitness-flux-analysis/scripts/forecast_similarity.py` against the production
HA1 clade map, so distances are identical across arms by construction.

## Results, matched subset, 180-day endpoint

Matched = identical region-weeks and forecast origins in both arms.

| arm | MLR error | vs naive | win rate |
|---|---|---|---|
| 2-year | 3.80 aa | +13.0% | 67.6% |
| 1-year | 3.89 aa | +9.0% | 59.8% |

Paired difference in MLR error (1-year minus 2-year), bootstrapped over forecast
origins:

| track | endpoint | difference | 95% CI |
|---|---|---|---|
| 365 d | 90 d | +0.074 aa | [−0.086, +0.268] |
| 365 d | 180 d | +0.090 aa | [−0.453, +0.558] |
| 365 d | 365 d | −0.217 aa | [−0.963, +0.600] |
| 180 d | 90 d | +0.054 aa | [−0.090, +0.228] |
| 180 d | 180 d | −0.032 aa | [−0.533, +0.408] |

The MLR-over-naive advantage gap (1-year minus 2-year) is likewise
indistinguishable: −3.95 pp [−15.4, +6.7] at the 180-day endpoint.

## Three things that would have produced a false positive

**1. Start-anchoring the windows.** A window's forecast origin is its *final* data
date, so pairing windows by equal label has the 1-year arm forecasting from a year
earlier than the 2-year arm — window length confounded with calendar epoch. The
1-year window labelled L must be paired with the 2-year window labelled L−1y, so
both end on the same day. This was got wrong first time round and is easy to get
wrong again.

**2. Scoring each arm on whatever it retained.** Unmatched, 1-year looked *better*
at 180 days (3.94 vs 4.04 aa); matched, it is worse (3.89 vs 3.80). The apparent
win came entirely from which region-weeks survived `location_min_seq` — the
1-year arm keeps only the data-rich ones, which are easier to forecast. Any future
comparison must match on `(region, forecast origin, truth origin, date, lead)`.
Note that matching on window *label* silently yields zero overlap, because the
labels differ by a year under end-anchoring.

**3. Treating region-weeks as independent.** There are only 18–21 independent
forecast origins; observations within one are strongly correlated. A per-observation
test over ~300 rows calls these differences significant. A bootstrap clustered on
forecast origin does not.

## Caveats on the result itself

The 1-year arm names a coarser clade set — 10.4 clades per window against 13.5,
missing 3.3 of the 2-year arm's clades on average, with more mass in the unscored
`other` bucket. Per `inclusion-thresholds/clade-analysis.md`, coarsening
mechanically lowers frequency error, and scoring in HA1 distance does not fix it
(that removes label exchangeability, not partition coarseness). So the 1-year arm
had a tailwind and still did not win, which strengthens the conclusion.

Do not read the label-based MAE column across arms: the 1-year arm looks much
worse there (8.6% vs 6.6%), but that is largely the divisor, since fewer named
clades means a smaller union to average over.

Three 1-year windows (2020-04, 2020-07, 2020-10) cannot be fit at all — no region
clears `location_min_seq: 500` over a 1-year COVID-trough span, which trips a hard
assert in `scripts/prepare-data.py`. Two more (2019-10, 2020-01) survive with only
two regions, where hierarchical pooling is nearly inert.

## To redo it

Generate 1-year window blocks end-anchored to their 2-year partners into a config
overlay (`extra_datasets`, `extra_analyses`, a `generation_time` entry for the new
virus token, and `clade_distance.map_source` pointing at the production map). The
dataset name must differ in the token *before the first underscore* — anything
like `h3n2_clades_1y_*` is silently absorbed into the `h3n2_clades` analysis by
`ff_io.seasonal_timepoints`, which globs `mlr-estimates/h3n2_clades_*`, and then
crashes `timepoint_to_numeric`. Wiring needed: extend `config["datasets"]` from
`extra_datasets` in the `Snakefile`, widen the `{analysis}` wildcard constraint in
`rules/fitness_flux_analysis.smk`, and make the clade-map inputs of
`rule forecast_similarity` resolve through `clade_distance.map_source`.

33 fits take about 15 minutes at 8 cores.
