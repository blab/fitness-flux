---
title: Variant frequency dynamics and short-term forecasting in SARS-CoV-2 and seasonal influenza
authors:
  - name: Trevor Bedford
    affiliations: [fredhutch, hhmi]
affiliations:
  - id: fredhutch
    name: Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA
  - id: hhmi
    name: Howard Hughes Medical Institute, Seattle, WA, USA
date: 2026-07-06
version:
  name: v1.0.1
  url: https://github.com/blab/fitness-flux/releases/tag/v1.0.1
links:
  - name: github.com
    url: https://github.com/blab/fitness-flux
  - name: bedford.io
    url: https://bedford.io/papers/bedford-fitness-flux
  - name: biorxiv.org
    url: https://www.biorxiv.org/content/10.64898/2026.07.05.736619
  - name: doi.org
    url: https://doi.org/10.64898/2026.07.05.736619
citation: "Bedford T. 2026. Variant frequency dynamics and short-term forecasting in SARS-CoV-2 and seasonal influenza. bioRxiv: 2026.07.05.736619."  
---

## Abstract

RNA virus populations are continually restructured by the emergence of new variants and their replacement of existing diversity, and tracking these variant frequencies is central to genomic surveillance.
Here we use multinomial logistic regression (MLR) to estimate the frequencies and relative growth rates of co-circulating variants of SARS-CoV-2 and the seasonal influenza lineages H3N2, H1N1pdm and B/Victoria, fitting the model in sliding windows across years of sequence data.
Beyond describing this turnover retrospectively, the fitted growth rates project variant frequencies forward, raising the question of how accurate those projections are as short-term forecasts.
Evaluated against retrospective truth over a six-month horizon, MLR forecasts of named-clade frequencies beat a naive persistence baseline for SARS-CoV-2 (mean absolute error 3.4% versus 5.0%) and for most seasonal influenza (H3N2 4.6% versus 5.8%; B/Victoria 6.2% versus 6.6%).
Frequency-based modeling thus offers both a compact description of ongoing variant turnover and a transparent baseline for anticipating it.

_The website [blab.github.io/fitness-flux/](https://blab.github.io/fitness-flux/) is the intended reading experience of this paper, providing responsive layout and interactive figures._

## Introduction

RNA viruses evolve rapidly and the selective pressures they face shift over the course of emergence to endemicity.
A newly established virus initially adapts to its new host by refining its capacity for within-host replication and between-host transmission.
Once a virus becomes endemic, adaptation is instead dominated by continual escape from accumulating population immunity, driving ongoing antigenic change and the recurrent emergence and replacement of variants [@kistler2023atlas].
This turnover is what genomic surveillance observes most directly: the rise and fall of variant frequencies through time.

Multinomial logistic regression (MLR) models the frequencies of co-circulating variants through time and infers a relative growth rate for each [@obermeyer2022analysis; @abousamra2024fitness].
Because it expresses this advantage as a difference in growth rate between variants, the measure maps directly onto the population-genetic notion of selective advantage, and these growth-rate differences correspond to differences in the time-varying effective reproduction number between co-circulating variants [@figgins2025frequency].
This frequency-based view of fitness [@volz2023fitness] has been widely applied to SARS-CoV-2 and provides a compact, interpretable description of variant dynamics.

Here we use MLR to trace the frequency dynamics of SARS-CoV-2 from the early pandemic in 2020 through 2025, alongside the three seasonal influenza lineages that circulate in humans, H3N2, H1N1pdm and B/Victoria [@bedford2014integrating].
Because the fitted growth rates also project variant frequencies forward, we then ask how well MLR forecasts short-term clade frequencies relative to a naive persistence baseline.

## Results and discussion

### Frequency dynamics

We follow population genetics first principles to compute the frequency through time of a haploid allele under deterministic selection in a fixed-generation Wright–Fisher population.
If an allele is at frequency $x$, then after a single generation with selective advantage $s$ the expected allele frequency will be 
$$x' = \frac{x \, (1+s)}{x \, (1+s) + (1-x)}.$$
Compounded over $t$ generations, the expectation from initial frequency $p$ follows
$$x(t) = \frac{p \, (1+s)^t}{p \, (1+s)^t + (1-p)}.$$
Generalizing this two-allele model to $n$ co-circulating variants, each with initial frequency $p_i$ and selective advantage $s_i$, variant $i$'s frequency is its relative abundance normalized by the sum across all variants,
$$x_i(t) = \frac{p_i \, (1+s_i)^t}{\sum_j p_j \, (1+s_j)^t}.$$
Moving from many discrete generations to continuous time, $(1+s_i)^t \approx \mathrm{exp}(t \, \mathrm{log}(1+s_i))$, so writing the growth rate $f_i = \mathrm{log}(1+s_i)$ gives the probability that a virus sampled at time $t$ is labeled as variant $i$
$$\mathrm{Pr}(X = i) = x_i(t) = \frac{p_i \, \mathrm{exp}(f_i \, t)}{\sum_j p_j \, \mathrm{exp}(f_j \, t) }.$$

This is the multinomial logistic regression (MLR) model, which has been widely used for modeling SARS-CoV-2 variant frequencies [@obermeyer2022analysis; @abousamra2024fitness].
The denominator normalizes the exponential growth/decay of individual variants so that overall frequency sums to 1.
The model has $2n$ parameters, with each variant $i$ having an initial frequency $p_i$ and a fixed growth rate $f_i$.
Because growth rates are necessarily relative, we fix an arbitrary "pivot" variant as a reference with growth rate $f=0$.
MLR growth rates are directly estimated in terms of calendar time with per-day or per-year values of $f$; we refer to this relative growth rate as the fitness of variant $i$.

We estimate frequencies of SARS-CoV-2 clades in 1-year sliding windows between Jan 2020 and Dec 2025 ([@fig:time-vs-frequency-sarscov2]).
In each window we collect clade sequence counts for viruses sampled from the USA and estimate per-variant frequencies and growth rates.
We use only the USA, the one country with sufficient temporal sequencing coverage over this period.
We collapse rare clades into a single "other" clade for MLR analysis to prevent noisy estimates from low sequence counts (see Methods).
The match between the empirical frequencies (dotted trajectories) and MLR frequencies (solid trajectories) indicates the model fits well despite having few parameters.

:::figure{#fig:time-vs-frequency-sarscov2 component=time-vs-frequency dataset=sarscov2_clades}
**Relative frequencies of SARS-CoV-2 clades through time.**
Points represent empirical frequencies of SARS-CoV-2 Nextstrain clades, while solid lines represent modeled frequencies from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
MLR frequency lines are drawn where there is sufficient sequence data to estimate empirical frequencies.
:::

For comparison purposes, we take a similar approach to estimating frequencies of seasonal influenza H3N2 ([@fig:time-vs-frequency-h3n2]).
Here we use 2-year sliding windows to account for slower frequency dynamics in seasonal influenza and still only use data from the USA.
The model fits are worse for H3N2 compared to SARS-CoV-2.
This is especially apparent at junctions between influenza seasons where stochastic seeding of a new season may result in a discontinuity of clade frequency compared to MLR expectation.
However, H3N2 fits remain sufficient to characterize clade frequency dynamics.

:::figure{#fig:time-vs-frequency-h3n2 component=time-vs-frequency dataset=h3n2_clades}
**Relative frequencies of H3N2 clades through time.**
Points represent empirical frequencies of H3N2 Nextstrain clades, while solid lines represent modeled frequencies from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
MLR frequency lines are drawn where there is sufficient sequence data to estimate empirical frequencies.
:::

### Seasonal influenza H1N1pdm and B/Victoria

We apply the identical frequency analysis to the two other seasonal influenza lineages that circulate in humans, H1N1pdm and B/Victoria, over the same 2016 to 2025 period, using USA data, 2-year sliding windows advanced every three months, a 2-week empirical-frequency smoothing window, and Nextclade HA subclades (`subclade_nextclade_ha`) as the variant classification.
Windows resolving fewer than two fittable USA clades are omitted, including those falling entirely within the 2020-2022 low-circulation trough when seasonal influenza all but vanished, so the series carry a gap across that period (roughly 2020 to late 2022 for H1N1pdm and 2020 to early 2023 for B/Victoria); H1N1pdm additionally lacks HA subclade resolution before ~2019, where its series begins.

Relative frequencies for both lineages show the same alternation of clade emergence and replacement seen in H3N2 ([@fig:time-vs-frequency-h1n1pdm]; [@fig:time-vs-frequency-vic]).

:::figure{#fig:time-vs-frequency-h1n1pdm component=time-vs-frequency dataset=h1n1pdm_clades}
**Relative frequencies of H1N1pdm clades through time.**
Points represent empirical frequencies of H1N1pdm Nextstrain clades, while solid lines represent modeled frequencies from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
MLR frequency lines are drawn where there is sufficient sequence data to estimate empirical frequencies.
:::

:::figure{#fig:time-vs-frequency-vic component=time-vs-frequency dataset=vic_clades}
**Relative frequencies of B/Victoria clades through time.**
Points represent empirical frequencies of B/Victoria Nextstrain clades, while solid lines represent modeled frequencies from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
MLR frequency lines are drawn where there is sufficient sequence data to estimate empirical frequencies.
:::

### Forecasting accuracy

Across all four lineages the MLR fits provide not only retrospective frequency estimates but forward projections of clade frequency, raising the question of how accurate those projections are as forecasts.
We adapt the real-time forecast-evaluation framework of Abousamra et al. [@abousamra2024fitness] to our sliding-window design ([@fig:forecast-accuracy]).
Treating each window's final date as the date of estimation, we forecast clade frequencies over the following six months and compare two models: MLR, which projects each named clade forward under its fitted logistic growth, and a naive model, which simply holds the final-date frequencies constant.
Because all clade windows now advance every three months, we take the window fit six months later (two quarters ahead) as the source of retrospective truth — its empirical smoothed frequency over the forecast period — and summarize accuracy as the mean absolute error (MAE) between predicted and observed frequencies, averaged across clades and across window pairs as a function of forecast lead time.
We score only the named clades, renormalized so they sum to one, and exclude the aggregated "other" bucket: it is a heterogeneous mix of rare and newly emerging clades that no fixed-growth model can anticipate, and retaining it otherwise dominates the long-lead error without reflecting predictive skill.

Over the three months preceding the date of estimation the two models coincide by construction — both are the in-window MLR fit — and track the retrospective frequencies to within 1-4% MAE.
They diverge only once projected past the date of estimation, and there the MLR projection beats naive persistence for three of the four lineages.
For SARS-CoV-2 the advantage is decisive (mean forecast MAE 3.4% for MLR versus 5.0% for naive), and it holds for H3N2 (4.6% versus 5.8%) and B/Victoria (6.2% versus 6.6%).
H1N1pdm is the exception: MLR (9.0%) roughly matches naive (8.4%) overall, winning at short lead but losing beyond a few months, consistent with its sparse, churny clade dynamics.
In every lineage both models rise with lead as clades that had not yet emerged at the date of estimation come to dominate the named set.
The absolute errors are smaller than the real-time figures reported by Abousamra et al. [@abousamra2024fitness], as expected, since we forecast from the fully backfilled retrospective fit rather than the sparse, delayed data available in real time, and our naive baseline is the MLR nowcast at the date of estimation rather than a raw recent-frequency average.
Overall, the fitted MLR growth advantages carry genuine short-term predictive signal for named-clade frequencies — pronounced for SARS-CoV-2 and present for most seasonal influenza.

:::figure{#fig:forecast-accuracy component=forecast-accuracy dataset=all}
**Forecasting accuracy of MLR versus a naive model across the four lineages.**
Mean absolute error between predicted and observed named-clade relative frequencies (the aggregated "other" bucket excluded and named clades renormalized to sum to one) as a function of forecast lead time, averaged across clades and across window pairs, for each virus lineage.
Negative leads (hindcast) fall within the fitting window, where the MLR and naive predictions coincide (grey); positive leads (forecast) project past each window's final date, where the MLR projection (blue) is compared to naive persistence (red).
The dashed line marks 5% error; faint lines show the individual window pairs, and $n$ is the number of pairs contributing to each panel.
:::

### Scoring forecasts by amino acid distance rather than clade label

The error above treats clade labels as exchangeable: a forecast that puts its mass on a clade's immediate parent is charged exactly as much as one that puts it on a distant clade.
For vaccine strain selection that is the wrong accounting, because a vaccine matched to a near relative of the clade that sweeps is a near hit while one matched to a distant clade is a miss.
The distinction is not hypothetical for H3N2: every one of the 40 windows contains at least one parent-descendant clade pair co-circulating in its named set.

To score forecasts in a way that tracks this, we represent each clade by the reconstructed HA1 amino-acid sequence at its most recent common ancestor and measure the distance between clades as the number of HA1 positions at which those sequences differ ([@fig:clade-distance]).
Across the 48 clades in the resulting map this spans 60 variable HA1 positions, with pairwise distances from 0 to 34 amino acids.
Clade labels are a lossy stand-in for that quantity: ancestor-descendant pairs differ by 3.8 amino acids on average and sibling pairs by 3.9, against 15.9 for pairs whose labels share no ancestry, but the three distributions overlap across the whole range — and J and J.3, which the label-based score treats as entirely distinct, have identical HA1 MRCA sequences.

We then score the *same* forecasts in the space of HA1 mutation profiles: each clade-frequency vector becomes, for every variable HA1 position, the fraction of the population carrying each amino acid, and the error is the distance between the forecast profile and the observed profile ([@fig:similarity-accuracy]).
Because a single global map covers every clade in the tree, the forecast and the truth no longer need to share a clade set — each is renormalized over its own window's named clades, and clades that emerged after the date of estimation enter the truth at their own coordinates rather than being dropped.
Error is then read in amino acids, and nothing about the forecasts themselves changes.

Scored this way MLR's advantage over persistence is larger than the clade-label score reports, but it is concentrated at intermediate lead rather than sustained to the end of the year.
Averaged over the whole one-year forecast window, MLR forecasts are off by 3.85 HA1 amino acids against 4.10 for persistence, a 6.3% advantage, while the same forecasts scored by clade label differ by 2.1% (6.25% against 6.39% mean absolute error).
Table 1 breaks that out at the origin and four forecast horizons.
The nowcast anchors the scale: at the origin both models reduce to the same in-window fit, off by 0.63 HA1 amino acids of fitting and sampling noise, and the MLR projection first separates from persistence at one month (10.9% closer by HA1 distance).
At three and six months MLR is clearly ahead by HA1 distance, with mean error 10.9% and 8.5% below persistence, while the clade-label score credits it with 7.1% and 2.1%.
By twelve months the advantage has gone under both scores: 2.0% by HA1 distance and 0.2% by clade label, which is to say the two models are indistinguishable.

The win rate — the fraction of region-weeks in which MLR lands strictly closer than persistence — does not follow the mean-error advantage.
At three months MLR's mean error is 10.9% below persistence by HA1 distance and yet it wins only 50.9% of region-weeks, so the early gain comes from a minority of windows in which MLR is far closer rather than from its being closer more often.
The win rate instead peaks at six months, at 62.6% by HA1 distance against 52.9% by clade label, and falls to 45.4% at twelve months against 32.6% by label.
The amino-acid score therefore does not rescue MLR at the longest lead.
What it shows is that MLR's six-month edge is real and substantially larger than the label score credits, and that by a year neither model holds a useful edge over the other.
Where MLR does win, it wins because when it is wrong about which clade will dominate it is wrong in the direction of a near relative, while persistence — anchored on the currently dominant clade — more often sits further away in HA1 space.
The label score is also increasingly unable to tell the two models apart at all: it returns exactly equal errors in 0.9% of region-weeks at three months, 3.4% at six and 24.6% at twelve, where HA1 distance never ties.

**Table 1. H3N2 forecast accuracy at the origin and at one, three, six and twelve months, scored by HA1 distance and by clade label.**
The same MLR and naive persistence forecasts, compared against the same observed frequencies, differing only in how the distance between a forecast and the truth is measured.
The clade-label rows are the mean absolute frequency error used elsewhere in this paper, taken over the union of the forecast's and the truth's clade sets since the two need not agree; that divisor is common to both models and so affects only the level, not the comparison.
MLR advantage is the reduction in mean error relative to naive; win rate is the fraction of region-weeks in which MLR's error is strictly smaller, with ties counted against MLR.
At the origin (nowcast) both models are the same in-window fit, so their errors coincide and only the error level is informative — a floor of 0.63 HA1 amino acids of fitting and sampling noise, against which the forecast rows should be read; the advantage and win-rate comparison is undefined there (dashes).
The 30- through 365-day blocks come from the one-year pairing (250 window pairs), covering leads of 25 to 30, 85 to 90, 175 to 180 and 360 to 365 days; the nowcast is taken from the six-month pairing, whose successor window covers the origin.
The clade-label score returns exact ties in 1.1%, 0.9%, 3.4% and 24.6% of region-weeks at 30, 90, 180 and 365 days respectively; the HA1 distance never ties.

| Score | MLR | Naive | MLR advantage | MLR win rate |
|---|---|---|---|---|
| *At the origin (nowcast)* | | | | |
| HA1 distance | 0.63 aa | 0.63 aa | — | — |
| Clade label (mean absolute error) | 1.43% | 1.43% | — | — |
| *At the 30-day endpoint* | | | | |
| HA1 distance | 1.06 aa | 1.19 aa | 10.9% | 55.4% |
| Clade label (mean absolute error) | 2.30% | 2.47% | 6.6% | 48.1% |
| *At the 90-day endpoint* | | | | |
| HA1 distance | 2.45 aa | 2.75 aa | 10.9% | 50.9% |
| Clade label (mean absolute error) | 4.30% | 4.63% | 7.1% | 50.6% |
| *At the 180-day endpoint* | | | | |
| HA1 distance | 4.04 aa | 4.42 aa | 8.5% | 62.6% |
| Clade label (mean absolute error) | 6.85% | 6.99% | 2.1% | 52.9% |
| *At the 365-day endpoint* | | | | |
| HA1 distance | 6.47 aa | 6.60 aa | 2.0% | 45.4% |
| Clade label (mean absolute error) | 9.00% | 9.02% | 0.2% | 32.6% |

The advantage is not uniform across regions: at one year it is positive in seven of ten regions, and negative in Africa, Europe and South Asia, the last of which contributes the fewest sequences.

:::figure{#fig:clade-distance component=clade-distance dataset=h3n2_clades}
**HA1 distances between H3N2 clades, and how poorly clade labels track them.**
Left, pairwise distance between the reconstructed HA1 amino-acid sequences at each clade's most recent common ancestor, for the 48 clades in the H3N2 similarity map, ordered by divergence from the root so that the tree's block structure is visible.
Right, the same 1128 clade pairs split by what their labels alone imply about relatedness, with the group mean marked.
The three distributions overlap across the full 0-34 amino acid range, so label identity is a weak proxy for HA1 similarity; a label-based forecast score charges the same penalty for every pair shown.
:::

:::figure{#fig:similarity-accuracy component=similarity-accuracy dataset=h3n2_clades}
**The same H3N2 forecasts scored by HA1 distance and by clade label.**
Left, forecast error in HA1 amino acids against lead time over a one-year horizon, for MLR (blue) and naive persistence (red); negative leads fall within the fitting window where the two coincide (grey), and faint lines show individual window pairs.
Right, the identical forecasts scored by clade label as mean absolute frequency error, where the two models nearly converge.
Both panels come from the same forecasts and the same observed frequencies; only the distance between them differs.
:::

The same construction applies to SARS-CoV-2 using spike in place of HA.
We build one global clade map from the Nextstrain all-time ncov tree, represent each Nextstrain clade by the reconstructed spike amino-acid sequence at its most recent common ancestor, and restrict the distance to the S1 subunit (positions 14 to 685) as the antigenic analog of HA1 ([@fig:clade-distance-sarscov2]).
Across the 56 clades in the map this spans 105 variable S1 positions, with pairwise distances from 0 to 59 amino acids.
Because Nextstrain clade labels such as 21K, 22B and 23A are not hierarchically nested, we read relatedness from the tree topology rather than the label string: ancestor-descendant clade pairs differ by 21.8 S1 amino acids on average and sibling pairs by 9.8, against 30.6 for unrelated pairs, but as for H3N2 the three distributions overlap across the whole range, and 17 clade pairs have identical S1 MRCA sequences.

Scored in spike S1 mutation-profile space, MLR again beats persistence by more than the clade-label score reports, and the two scores diverge more sharply than for H3N2 ([@fig:similarity-accuracy-sarscov2]).
Averaged over the one-year forecast window MLR forecasts are off by 10.0 S1 amino acids against 10.8 for persistence, a 7.6% advantage, while the same forecasts scored by clade label differ by 3.1% (6.81% against 7.02% mean absolute error).
Table 2 breaks this out at the origin and four forecast horizons.
The nowcast again anchors the scale — both models reduce to the same fit, off by 0.85 S1 amino acids — and at one month MLR is far ahead by both scores (37.5% closer by S1 distance and 25.7% by clade label, winning about 69% of region-weeks either way), the one horizon at which the label score still agrees.
At three, six and twelve months MLR remains ahead by S1 distance, with mean error 16.6%, 7.6% and 1.4% below persistence, and it lands strictly closer in a majority of region-weeks throughout (65.7%, 57.0% and 63.5%).
The clade-label score, by contrast, credits MLR with 11.9% at three months, exactly 0.0% at six and −0.6% at twelve, and is unable to tell the two models apart in a rising fraction of region-weeks — tied in 10.0%, 26.9% and 61.5% respectively, where the S1 distance never ties.
Read by clade label the two models look indistinguishable by six months and persistence marginally better by a year; read in spike space MLR keeps a small but consistent edge throughout — the same qualitative correction as for H3N2, and starker here because SARS-CoV-2's clade turnover moves more amino acids per step.

**Table 2. SARS-CoV-2 forecast accuracy at the origin and at one, three, six and twelve months, scored by spike S1 distance and by clade label.**
The same MLR and naive persistence forecasts, compared against the same observed frequencies, differing only in how the distance between a forecast and the truth is measured, exactly as in Table 1.
At the origin (nowcast) both models are the same in-window fit and their errors coincide — a floor of 0.85 spike S1 amino acids — so the advantage and win-rate comparison is undefined there (dashes).
The 30- through 365-day blocks come from the one-year pairing (80 window pairs across six regions), covering leads of 25 to 30, 85 to 90, 175 to 180 and 360 to 365 days; the nowcast is taken from the six-month pairing, whose successor window covers the origin (the one-year pairing's successor opens only after it, so it carries no nowcast).
The clade-label score returns exact ties in 4.5%, 10.0%, 26.9% and 61.5% of region-weeks at 30, 90, 180 and 365 days respectively; the spike S1 distance never ties.

| Score | MLR | Naive | MLR advantage | MLR win rate |
|---|---|---|---|---|
| *At the origin (nowcast)* | | | | |
| Spike S1 distance | 0.85 aa | 0.85 aa | — | — |
| Clade label (mean absolute error) | 1.30% | 1.30% | — | — |
| *At the 30-day endpoint* | | | | |
| Spike S1 distance | 1.26 aa | 2.02 aa | 37.5% | 68.7% |
| Clade label (mean absolute error) | 1.66% | 2.24% | 25.7% | 69.4% |
| *At the 90-day endpoint* | | | | |
| Spike S1 distance | 5.56 aa | 6.67 aa | 16.6% | 65.7% |
| Clade label (mean absolute error) | 5.23% | 5.94% | 11.9% | 50.5% |
| *At the 180-day endpoint* | | | | |
| Spike S1 distance | 10.89 aa | 11.79 aa | 7.6% | 57.0% |
| Clade label (mean absolute error) | 8.20% | 8.20% | 0.0% | 25.6% |
| *At the 365-day endpoint* | | | | |
| Spike S1 distance | 19.60 aa | 19.88 aa | 1.4% | 63.5% |
| Clade label (mean absolute error) | 9.23% | 9.17% | −0.6% | 7.1% |

:::figure{#fig:clade-distance-sarscov2 component=clade-distance dataset=sarscov2_clades}
**Spike S1 distances between SARS-CoV-2 clades, and how poorly clade labels track them.**
Left, pairwise distance between the reconstructed spike S1 amino-acid sequences at each Nextstrain clade's most recent common ancestor, for the 56 clades in the SARS-CoV-2 similarity map, ordered by divergence from the root so that the tree's block structure is visible.
Right, the same 1540 clade pairs split by tree-topology relatedness (ancestor/descendant, sibling, unrelated), with the group mean marked.
The three distributions overlap across the full 0-59 amino acid range, so relatedness is a weak proxy for S1 similarity; a label-based forecast score charges the same penalty for every pair shown.
:::

:::figure{#fig:similarity-accuracy-sarscov2 component=similarity-accuracy dataset=sarscov2_clades}
**The same SARS-CoV-2 forecasts scored by spike S1 distance and by clade label.**
Left, forecast error in spike S1 amino acids against lead time over a one-year horizon, for MLR (blue) and naive persistence (red); negative leads fall within the fitting window where the two coincide (grey), and faint lines show individual window pairs.
Right, the identical forecasts scored by clade label as mean absolute frequency error, where the two models converge and, past six months, are tied in most region-weeks.
Both panels come from the same forecasts and the same observed frequencies; only the distance between them differs.
:::

## Conclusions

Genomic surveillance rests on tracking which variants are rising and falling in frequency.
Modeling these frequencies with multinomial logistic regression provides both a compact description of variant turnover and, through the fitted per-variant growth rates, forward projections of frequency.
We find these projections carry genuine short-term predictive signal: across SARS-CoV-2 and seasonal influenza, MLR forecasts of named-clade frequencies beat naive persistence over a six-month horizon for all but the sparsest lineage, and most decisively for SARS-CoV-2.
Because the growth rates carry an absolute scale, they also connect back to differences in the effective reproduction number between co-circulating variants and to epidemiological impact [@figgins2025frequency].

The generality of the approach rests on a single requirement: a way to bin genetic diversity into discrete, comparable variants.
For SARS-CoV-2 and influenza this comes off the shelf, with Nextstrain clades supporting both the frequency and the forecast analysis.
Pathogens without an established nomenclature could be analyzed via automated methods that partition a tree into lineages [@mcbroome2024framework; @lefrancq2025learning].
These frequency-based forecasts provide a transparent baseline for anticipating variant success, against which more elaborate models can be measured.

## Methods

### Sequence data

For SARS-CoV-2, we use curated "open" data from Nextstrain [@hadfield2018nextstrain] that draws from NCBI GenBank.
For influenza H3N2, we use data from GISAID [@shu2017gisaid].
In each case, the raw sequences are processed with Nextclade [@aksamentov2021nextclade] to assign Nextstrain clade to SARS-CoV-2 sequences and to assign subclade [@neher2026nomenclature] to influenza sequences.
We filter out sequences with Nextclade overall QC status of "bad".
We additionally filter to sequences collected from the USA.
This leaves 3,588,802 total sequences for SARS-CoV-2 sampled between 2020 and 2025 and 44,456 total sequences for H3N2 sampled between 2016 and 2025.

### Multinomial logistic regression

We conducted multinomial logistic regression (MLR) using the evofr package ([github.com/blab/evofr](https://github.com/blab/evofr)) on sliding windows advanced every three months: 1-year windows for SARS-CoV-2 (23 windows) and 2-year windows for the seasonal influenza lineages (40 for H3N2, 28 for H1N1pdm, 35 for B/Victoria), keeping only windows that resolve at least two fittable clades.
For each window we treat each clade as a distinct variant, collapsing rare clades together into a single "other" category before fitting.
For all clade datasets, a clade is modeled separately only if it reaches at least 50 sequences and a mean frequency of at least 0.1% within the window, while clades below either threshold are merged into "other".
This leaves between 7 and 18 clades per window (median 14) for SARS-CoV-2, between 3 and 16 (median 9) for H3N2, and between 2 and 16 (median 6) and 2 and 8 (median 3) for the sparser H1N1pdm and B/Victoria lineages.
We set no explicit per-window pivot; evofr references the "other" bucket.
The "other" bucket is fit as a variant but excluded from all downstream analyses of clade frequency (the per-season frequencies and the forecast evaluation), which are computed over the named clades renormalized to sum to one.

We use generation time $\tau$ of 5.0 days for pre-Omicron SARS-CoV-2 following [@ferretti2020quantifying; @ganyani2020estimating; @hart2022generation], generation time of 3.2 days for post-Omicron SARS-CoV-2 following [@park2023inferring; @chan2026estimating] and generation time of 3.2 days for seasonal influenza H3N2 following [@cowling2009estimation; @vink2014serial; @chan2025estimating].
At first order, the fitness $f_i = \tau \, f_i^{\mathrm{day}}$ is the log ratio of reproduction numbers between variant $i$ and the pivot, $\mathrm{log}(R_i / R_{\mathrm{pivot}})$.
Only the mean generation interval $\tau$ enters at this order.
The shape of the generation-interval distribution enters at second order, informed by the interval's variance and absolute per-day growth rates rather than $f_i^{\mathrm{day}}$.
Because the frequency data identify only this relative rate, the second order correction cannot be estimated from frequency dynamics alone.
The per-generation multiplicative fitness $\mathrm{exp}(f_i) = 1 + s_i$ therefore equals the reproduction-number ratio $R_i / R_{\mathrm{pivot}}$ under a fixed (point) generation interval and is an upper bound under a realistically dispersed interval [@wallinga2006generation], for which the fuller frequency-to-reproduction-number conversion is given in Figgins and Bedford [@figgins2025frequency].

### Forecast accuracy

To evaluate forecasting accuracy ([@fig:forecast-accuracy]) we pair each window with the window fit two quarters (~180 days) later, which supplies retrospective truth over a six-month forecast horizon; pairing keys off the actual final data dates (not the nominal window bounds, since sparse influenza data can end early), so a window whose successor is missing, for instance across the low-circulation gap, is skipped.
Writing $T$ for the earlier window's final date and $S$ for its named clades (the "other" bucket excluded and $S$ renormalized to sum to one), we evaluate, at each date $t$ in $[T-90, T+180]$ days (a three-month hindcast and six-month forecast), three quantities over $S$.
The MLR forecast is $\hat{x}^{\mathrm{MLR}}_i(t) = \mathrm{softmax}_i[\log x_i(T) + f_i^{\mathrm{day}} (t - T)]$ for $t > T$, and the in-window modeled frequency for $t \le T$, where $f_i^{\mathrm{day}} = f_i / \tau_i$ is clade $i$'s per-day logistic growth rate recovered from its fitted (within-window) growth advantage.
The naive forecast is $\hat{x}^{\mathrm{naive}}_i(t) = x_i(T)$ held constant for $t > T$, and equal to the MLR fit for $t \le T$, so the two coincide over the hindcast and differ only in how they project forward.
The truth $x_i(t)$ is the empirical smoothed frequency from the later window, restricted to $S$ and renormalized over $S$ (mass in the later window's "other" bucket and in clades that emerged after $T$ is dropped).
The absolute error at each date is the mean over clades, $\frac{1}{|S|} \sum_{i \in S} \lvert x_i(t) - \hat{x}_i(t) \rvert$ [@abousamra2024fitness], which we average within weekly lead-time bins across each lineage's contributing window pairs (21 for SARS-CoV-2, 32 for H3N2, 19 for H1N1pdm and 22 for B/Victoria) to give mean absolute error as a function of lead time.

### Similarity-aware forecast scoring

To score forecasts by amino acid distance rather than clade identity ([@fig:clade-distance], [@fig:similarity-accuracy]) we build one global clade map for H3N2 from the Nextstrain seasonal-flu HA tree.
Walking the tree from the root and accumulating HA1 amino-acid mutations along each path gives the HA1 sequence at every internal node; a clade's representative sequence is the one at the shallowest node carrying that clade label, i.e. its most recent common ancestor.
We use the MRCA sequence rather than the consensus of the sequences assigned to a clade because the consensus is a mixture over whichever descendants happen to fall under that label and therefore shifts with the collapse threshold, while the MRCA sequence depends only on the tree.
The distance $d(k, l)$ between two clades is the number of HA1 positions at which their MRCA sequences differ.
One clade named in the analysis windows, B.1.2, has no node of its own in the tree and inherits its nearest reconstructable ancestor's sequence (B.1); `unassigned` has no MRCA and, unlike in the label-based score, is excluded here and the remaining frequencies renormalized (it averages 0.35% of H3N2 sequences).
For SARS-CoV-2 the identical procedure is applied to the Nextstrain all-time ncov tree with spike in place of HA, restricting the distance to the S1 subunit (positions 14 to 685); the tree's `clade_membership` labels are matched to the Nextstrain clade the frequency analysis names by their leading token (so "23A (XBB.1.5)" matches "23A"), and the `WT` bucket, which merges the basal 19A and 19B clades and carries no node of its own, takes the Wuhan-Hu-1 root sequence.

Each globally variable HA1 position contributes one column per amino acid observed at it across all clades, giving a clade $\times$ column indicator matrix $M$ with a 1 where clade $k$ carries column $c$'s amino acid.
One column per allele, rather than per non-reference allele, is what keeps the scale consistent: two clades differing at a single position differ in exactly two columns whether that position is biallelic or not.
Writing $A$ for that matrix scaled by $1/2$, a frequency vector $p$ over clades maps to the population mutation profile $A^{\mathsf{T}} p$, and the error between a forecast $\hat{x}(t)$ and the truth $x(t)$ is

$$ E(t) = \lVert A^{\mathsf{T}} \hat{x}(t) - A^{\mathsf{T}} x(t) \rVert_1 . $$

The scaling makes this read directly in amino acids: two pure clades sit $d(k, l)$ apart.
The map is deliberately not injective — clades whose MRCA sequences are identical (H3N2 has one such pair, J and J.3) occupy the same point and the score cannot separate them, which is the intended behaviour for a measure of amino acid difference rather than of nomenclature.
Predictions are unchanged from the label-based score — the same MLR projection and the same persistence baseline, from the same fits — so the comparison isolates the scoring convention.
Because $M$ is global, the forecast is renormalized over the earlier window's named clades and the truth over the later window's, with no remapping between them; the label-based score by contrast restricts truth to the earlier window's clade set, which drops clades that emerged after $T$ entirely.
We score both the existing six-month pairing and a one-year pairing (each window paired with the one four quarters ahead), over $[T-90, T+h]$ days for horizon $h$.

The label-based comparator we report alongside is the mean absolute frequency error of the Forecast accuracy section, $\frac{1}{|U|} \sum_i \lvert x_i(t) - \hat{x}_i(t) \rvert$, where $U$ is the union of the forecast's and the truth's clade sets; the union stands in for the single clade set $S$ available there, because here the two sides need not agree.
That divisor is common to both models at a given date, so it sets the level of the error but cannot change which model is closer.
To verify the transform we add a uniform separation $\epsilon$ between every pair of distinct clades, via one identity column per clade, and confirm that as $\epsilon$ grows the normalized error converges to $\frac{1}{2}\sum_i \lvert x_i(t) - \hat{x}_i(t) \rvert$, the total-variation distance between the two frequency vectors and thus $\lvert U \rvert / 2$ times that mean absolute error — clade identity being all that survives once every pair is equidistant.
That check runs as a workflow step ahead of any scoring, together with the requirements that two pure clades sit exactly $d(k,l)$ apart and that moving forecast mass onto a more distant clade always raises the score.
Reported results use $\epsilon = 0$ throughout.

### Reproducibility

The workflow to reproduce data, analyses and manuscript is available at [github.com/blab/fitness-flux](https://github.com/blab/fitness-flux).
SARS-CoV-2 sequence data were sourced from publicly available open data in NCBI GenBank and so this side of the workflow is easily runnable.
However, influenza H3N2 sequence data was collected from GISAID and these data are not allowed to be reshared and consequently the influenza side of the workflow is unfortunately not reproducible without separately providing GISAID data.

## Acknowledgments

We thank Jumpei Ito, Michael Lässig, Erick Matsen, Richard Neher, Cécile Viboud and members of the Bedford Lab for helpful feedback.
SARS-CoV-2 analyses are based on open data in GenBank.
We gratefully acknowledge the researchers and data contributors who collected the specimens, generated and deposited the raw sequence data and metadata into NCBI GenBank.
Influenza analyses are based on GISAID data.
We gratefully acknowledge all data contributors, i.e., the Authors and their Originating laboratories
responsible for obtaining the specimens, and their Submitting laboratories for generating the
genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is
based. 
TB was funded as a Howard Hughes Medical Institute Investigator.
This research was supported in part by grant NSF PHY-2309135, the Gordon and Betty Moore Foundation grant no. 2919.02, and the Chan Zuckerberg Initiative DAF grant to the Kavli Institute for Theoretical Physics (KITP).
