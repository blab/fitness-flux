# H3N2 `location_min_seq` ladder, scored by HA1 distance

**Verdict: keep `location_min_seq: 500`.** Raising it does make the MLR-over-naive advantage
go up, significantly so at 1000 and 2000 — but the gain comes from persistence getting *worse*,
not from MLR getting better. MLR's own error is flat across the ladder, and the win rate is
flat too. Buying a wider margin by weakening the comparator, at the cost of half the
region-weeks and six of ten regions, is not a trade worth making.

The production config was never changed; the arms ran from a config overlay and have been
removed. This note is the whole artifact.

## What was run

Four arms over the production 2-year H3N2 window grid, identical in every respect except the
threshold: same dates, same `clade_min_count`/`clade_min_freq`, same hierarchical model and
`pool_scale: 0.25`, all scored against the production HA1 clade map so distances are identical
by construction. 116 new fits; the 500 arm is production and was not refit.

| `location_min_seq` | region-windows | regions/window | windows |
|---|---|---|---|
| 500 (current) | 354 | 8.85 | 40 |
| 1000 | 271 | 6.78 | 40 |
| 2000 | 181 | 4.53 | 40 |
| 4000 | 105 | 2.92 | 36 |

At 4000, `2019-04`, `2019-07`, `2019-10` and `2020-01` retain no region at all and trip the
hard assert in `scripts/prepare-data.py`; they were omitted from that arm. Region loss is
uneven — South Asia is gone by 2000, as are China and Japan Korea; only Europe and North
America reach 4000.

## Results at the 180-day endpoint

Each arm scored as configured, on whatever region-weeks it retains (365-day pairing track):

| `location_min_seq` | HA1 MLR | naive | advantage | win rate |
|---|---|---|---|---|
| 500 | 4.044 aa | 4.420 | 8.5% | 62.6% |
| 1000 | 3.954 | 4.480 | 11.7% | 61.0% |
| **2000** | 3.987 | 4.570 | **12.8%** | 62.4% |
| 4000 | 4.014 | 4.358 | 7.9% | 62.2% |

An inverted U peaking at 2000, reproduced in the 180-day pairing track (8.5 → 14.0 → 15.5 →
11.8%). Bootstrapped paired on forecast origin, the advantage gap over the 500 baseline is
real at the middle rungs:

| arm | 365 d track | 180 d track |
|---|---|---|
| 1000 vs 500 | +3.03 pp [+0.28, +6.45] | +3.25 pp [+0.48, +6.73] |
| 2000 vs 500 | +4.53 pp [+1.35, +8.79] | +4.77 pp [+1.64, +8.93] |
| 4000 vs 500 | +2.30 pp [−4.19, +11.81] | +3.90 pp [−3.12, +13.54] |

## Why that is not a reason to swap

Decomposing the ratio on shared forecast origins, 180-day endpoint:

| | MLR error change vs 500 | naive error change vs 500 |
|---|---|---|
| 1000, 365 d track | −0.068 aa [−0.215, +0.081] | +0.042 [−0.037, +0.174] |
| 2000, 365 d track | −0.003 aa [−0.314, +0.331] | +0.202 [−0.123, +0.620] |
| 1000, 180 d track | −0.081 aa [−0.221, +0.053] | +0.031 [−0.038, +0.147] |
| 2000, 180 d track | +0.034 aa [−0.184, +0.301] | **+0.257 [+0.019, +0.601]** |

**MLR's error does not move.** The only single-model change that clears zero anywhere is naive
getting worse. The advantage is a ratio, and it improves because the denominator degrades.

The win rate agrees and is the blunter statement: 62.6 / 61.0 / 62.4 / 62.2% across the ladder.
MLR does not beat persistence more *often* at a higher threshold — only by a wider margin on a
baseline that has become harder.

## The selection effect runs the opposite way to expectation

Going in, the worry was that raising the threshold would retain easier cells and flatter every
arm above the baseline. It does the reverse: **naive error rises with the threshold**, so the
retained cells are *harder* for persistence, not easier.

The likely reason is that the large regions (Europe, North America) have faster clade turnover,
while the sparse regions are comparatively static — and persistence does well precisely where
little changes. Dropping small regions removes the cells where doing nothing was a good
strategy. That is worth knowing on its own, and it is why the advantage metric misleads here
even though it is the better-controlled one.

## Things that had to be got right

**Pair the bootstrap on forecast origin.** Resampling the arms independently gave intervals of
roughly ±20 pp — wide enough to hide everything. All arms share the same window grid, so the
origins are a common set and the comparison can be paired, which is what made the effect
visible at all. There are only 16–20 shared origins at a given endpoint; per-observation
intervals over a few hundred correlated region-weeks would be several times too narrow.

**Lead with the advantage, but decompose it.** Absolute error is confounded by which cells each
arm retains; the advantage is not, since both models are scored on identical cells within an
arm. But "not confounded by selection" is not "measures skill" — a ratio can improve from
either end, and here it improved from the wrong one.

**The clade partition is invariant to this lever, so it is not a confound.**
`scripts/prepare-data.py` picks the clade keep-set (lines 189-218) over the location-*unfiltered*
frame and applies `locations_to_include` only at line 246. Confirmed empirically: mean named
clades were 13.53 / 13.50 / 13.40 / 13.69 across the four arms. This differs from the
window-length comparison, where coarsening was real — see
`notes/2026-09-23-window-length-comparison.md`.

**Forecast origins can shift between arms.** `metadata.dates` comes from the observed data, not
the config range, so dropping the region that supplied a window's last sequence moves
`dates[-1]` — which is the forecast origin, the zero point for every lead, and the
window-pairing key. It moved for 0 / 4 / 7 / 6 windows of 40. The analysis is restricted to
shared origins for this reason.

## Adjacent lever, untested

`--location-min-seq-days` restricts the counting window used to admit a region to the last N
days, and is set for **no dataset**. With 2-year windows, a region that sequenced heavily in
year one and then stopped still qualifies on 500 while contributing nothing near the forecast
origin. That is a better-targeted lever than raising the count — it would drop stale regions
rather than small ones — and is the natural next thing to try.

## To redo it

Generate an overlay (`extra_datasets`, `extra_analyses`, a `generation_time` entry per new
virus token, and `clade_distance.map_source` pointing every arm at the production map) with one
block per (arm, window) overriding only `location_min_seq`, and omitting windows that retain no
region. The dataset name must differ in the token *before the first underscore* — anything like
`h3n2_clades_loc2000_*` is silently absorbed into the production analysis by
`ff_io.seasonal_timepoints`, which globs `mlr-estimates/h3n2_clades_*`. Wiring: extend
`config["datasets"]` from `extra_datasets` in the `Snakefile`, widen the `{analysis}` wildcard
constraint in `rules/fitness_flux_analysis.smk`, and resolve `rule forecast_similarity`'s map
inputs through `clade_distance.map_source`.

116 fits take about 25 minutes at 8 cores.
