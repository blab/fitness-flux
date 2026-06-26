---
title: Fitness flux in SARS-CoV-2 and influenza H3N2
authors:
  - name: Trevor Bedford
    affiliations: [fredhutch, hhmi]
affiliations:
  - id: fredhutch
    name: Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA
  - id: hhmi
    name: Howard Hughes Medical Institute, Seattle, WA, USA
date: 2026-06-21
---

## Abstract

The tempo of viral adaptation is usually read indirectly from the composition of mutations, through measures such as dN/dS.
Here we measure it directly from the dynamics of variant frequencies, where we use multinomial logistic regression to estimate a fitness for each co-circulating variant.
We aggregate these estimates to derive the rate of change of mean population fitness, referred to as fitness flux.
Tracing SARS-CoV-2 from 2020 to 2026 and comparing against seasonal influenza A/H3N2, we find that SARS-CoV-2 adapted rapidly with a 6.7-fold increase in fitness from 2020 to 2023, but then slowing with a 2.2-fold increase in fitness from 2023 to 2026.
Influenza H3N2 sustains a slower, steadier pace roughly threefold below recent SARS-CoV-2.
In both, the rate of fitness gain closely tracks the variance in fitness, matching the 1:1 expectation of Fisher's fundamental theorem.
Phylogenetic contrasts between parent and child lineages localize most fitness gain to spike, and within spike to the receptor-binding domain, where a simple count of spike S1 substitutions predicts lineage fitness about as well as deep-learning escape and protein-language-model scores.
Measuring fitness directly thus offers a transparent, frequency-based alternative to mutational proxies for tracking and anticipating viral adaptation.

## Introduction

RNA viruses evolve rapidly and the selective pressures they face shift over the course of emergence to endemicity.
A newly established virus initially adapts to its new host by refining its capacity for within-host replication and between-host transmission.
Once a virus becomes endemic, adaptation is instead dominated by continual escape from accumulating population immunity, driving ongoing antigenic change, where some viruses are better able than others to sustain ongoing adaptive evolution [@kistler2023atlas].
Adaptation of either kind leaves a signature in amino-acid replacing nonsynonymous vs silent synonymous substitutions, with methods ranging from simple comparisons of nonsynonymous to synonymous substitution rates (dN/dS) to McDonald–Kreitman-style approaches [@mcdonald1991adaptive] that weigh mutations fixed along a virus's successful trunk lineage against those lost on unsuccessful side branches [@wolf2006long].
Such approaches have revealed rapid, continued adaptation in the SARS-CoV-2 spike S1 subunit [@kistler2022rapid; @markov2023evolution].

A complementary class of methods estimates fitness directly from the dynamics of variant frequencies rather than from the composition of mutations [@volz2023fitness].
Multinomial logistic regression (MLR) models the frequencies of co-circulating variants through time and infers a relative growth rate, or fitness, for each [@obermeyer2022analysis; @abousamra2024fitness].
Because it expresses fitness as a difference in growth rate between variants, this measure maps directly onto the population-genetic notion of selective advantage.
These growth-rate differences correspond to differences in the time-varying effective reproduction number between co-circulating variants [@figgins2025frequency].
Aggregating these per-variant fitnesses into the rate of change of mean population fitness yields the population's fitness flux [@mustonen2010fitness], a direct, frequency-based alternative to dN/dS for quantifying the tempo of adaptation.

Here we use this frequency-based view of fitness to trace how SARS-CoV-2 has adapted from the early pandemic in 2020 through to 2026, spanning the transition from initial host adaptation to sustained evolution for antigenic novelty.
We place the rate of SARS-CoV-2 fitness change in context by comparing it against seasonal influenza A/H3N2 which exhibits canonically rapid and continuous adaptation [@bedford2014integrating].
Finally, we relate the inferred changes in fitness to molecular predictors, most directly the accumulation of spike mutations, to identify the substitutions that drive fitness gain.

## Results and discussion

### Frequency dynamics

We follow population genetics first principles to compute the frequency through time of a haploid allele under selection.
If an allele is at frequency $x$, then after a single generation with selective advantage $s$ the expected allele frequency will be 
$$x' = \frac{x \, (1+s)}{x \, (1+s) + (1-x)}.$$
Compounded over $t$ generations, the expectation from initial frequency $p$ follows
$$x(t) = \frac{p \, (1+s)^t}{p \, (1+s)^t + (1-p)}.$$
Generalizing this two-allele model to $n$ co-circulating variants, each with initial frequency $p_i$ and selective advantage $s_i$, variant $i$'s frequency is its relative abundance normalized by the sum across all variants,
$$x_i(t) = \frac{p_i \, (1+s_i)^t}{\sum_j p_j \, (1+s_j)^t}.$$
Moving from many discrete generations to continuous time, $(1+s_i)^t \approx \mathrm{exp}(t \, \mathrm{log}(1+s_i))$, so writing the growth rate $f_i = \mathrm{log}(1+s_i)$ gives the probability that a virus sampled at time $t$ is labeled as variant $i$
$$\mathrm{Pr}(X = i) = x_i(t) = \frac{p_i \, \mathrm{exp}(f_i \, t)}{\sum_j p_j \, \mathrm{exp}(f_j \, t) }.$$

This is the multinomial logistic regression (MLR) model, which has seen significant previous use for modeling SARS-CoV-2 variant frequencies [@obermeyer2022analysis; @abousamra2024fitness].
The denominator normalizes the exponential growth/decay of individual variants so that overall frequency sums to 1.
The model has $2n$ parameters, with each variant $i$ having an initial frequency $p_i$ and a fixed growth rate $f_i$.
Because growth rates are necessarily relative, we fix an arbitrary "pivot" variant as a reference with growth rate $f=0$.
MLR growth rates are directly estimated in terms of calendar time with per-day or per-year values of $f$.
To express fitness in per-generation units we multiply each per-day rate by the generation time $\tau$ measured in days, giving $f_i = \tau \, f_i^{\mathrm{day}}$, the change in log frequency accrued over a single generation.
We assume $\tau = 3.2$ days for both SARS-CoV-2 and influenza H3N2 [@song2022serial].
Throughout, we refer to this per-generation growth rate $f_i = \mathrm{log}(1+s_i)$ as the fitness of variant $i$.
Because $f_i$ is defined on a log scale, mean fitness, fitness variance, fitness flux and changes in fitness between lineages are all likewise computed on this log scale.

We estimate frequencies and fitnesses of SARS-CoV-2 clades in 1-year sliding windows between Jan 2020 and Jan 2026 ([@fig:time-vs-frequency-sarscov2]).
In each window we collect clade sequence counts for viruses sampled from the USA and estimate per-variant frequencies and fitnesses.
We use just the USA in this analysis because the USA is the only country with sufficient temporal data with good sequencing coverage from 2020 through 2025.
We collapse rare clades into a single "other" clade for MLR analysis to prevent noisy estimates from low sequence counts (see Methods).
The match between the empirical frequencies (dotted trajectories) and MLR frequencies (solid trajectories) indicates the model fits well despite having few parameters.

:::figure{#fig:time-vs-frequency-sarscov2 component=time-vs-frequency dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_frequency.png}
**Relative frequencies of SARS-CoV-2 clades through time.**
Points represent empirical frequencies of SARS-CoV-2 Nextstrain clades, while solid lines represent modeled frequencies from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

For comparison purposes, we take a similar approach to estimating frequencies and fitnesses of seasonal influenza H3N2 ([@fig:time-vs-frequency-h3n2]).
Here we use 2-year sliding windows to account for slower frequency dynamics in seasonal influenza and still only use data from the USA.
The model fits are worse for H3N2 compared to SARS-CoV-2.
This is especially apparent at junctions between influenza seasons where stochastic seeding of new season may result in a discontinuity of clade frequency compared to MLR expectation.
However, we believe that H3N3 model fits are sufficient to correctly estimate the magnitude of fitness effects.

:::figure{#fig:time-vs-frequency-h3n2 component=time-vs-frequency dataset=h3n2_clades static=figures/h3n2_clades_time_vs_frequency.png}
**Relative frequencies of H3N2 clades through time.**
Points represent empirical frequencies of H3N2 Nextstrain clades, while solid lines represent modeled frequencies from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

We "scaffold" MLR fitness estimates across windows to arrive at a single fitness estimate per variant.
Here, each window only measures fitness differences, with its own arbitrary zero.
We solve for the single set of clade fitnesses and per-window offsets that best fit every window at once, weighting by abundance (see Methods).

### Fitness flux

With variant frequency $x_i(t)$ and constant variant fitness $f_i$, we describe the mean population fitness as a standard weighted sum $\bar{f}(t) = \sum_i x_i(t) \, f_i$.
The fitness flux [@mustonen2010fitness] of the population is then the rate of change of population fitness at a given time $\phi(t) = \Delta \bar{f}(t) / \Delta t$.
Integrating this rate gives the cumulative fitness flux
$$\Phi(t) = \int_{t_0}^{t} \phi(t') \, dt' = \bar{f}(t) - \bar{f}(t_0),$$
the total adaptive change accumulated along the population's trajectory.
Because variant fitnesses are estimated only relative to a pivot, an individual variant's scaffolded fitness is meaningful as a difference from a baseline rather than as an absolute value.
Chaining these locally-measured advantages across overlapping windows places variant $i$ at a cumulative fitness flux $\Phi_i = f_i - f_0$ relative to the founding variant, and the population sits at the frequency-weighted average $\Phi(t) = \sum_i x_i(t) \, \Phi_i$.

We find that SARS-CoV-2 accumulates fitness flux rapidly, over the course of 2020 through 2025, doubling in fitness every 1.5 years on average ([@fig:time-vs-fitness-sarscov2]).
After initial spread of D614G [@korber2020tracking] in 2020, we observe a lull, followed by rapid growth in fitness in 2021 and 2022 with initial VOCs, Omicron and initial Omicron sub-lineages [@roemer2023sars], and then a slower, more steady pace since 2024.
There is a mix of large jumps in fitness (in particular to BA.1, but also more recently to JN.1) and smaller, more gradual step change.

:::figure{#fig:time-vs-fitness-sarscov2 component=time-vs-fitness dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_fitness.png}
**Cumulative SARS-CoV-2 fitness flux.**
Empirical frequencies of SARS-CoV-2 clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

Seasonal influenza H3N2 shows a fundamentally similar pattern of emergence of new clades and their replacement of existing diversity.
However, H3N2 dynamics play out on a slower timescale ([@fig:time-vs-fitness-h3n2]).
Rather than doubling in fitness every 1.5 years, H3N2 shows a doubling period of 9.7 years.
Greater coexistence of multiple co-circulating clades is also apparent relative to SARS-CoV-2.

:::figure{#fig:time-vs-fitness-h3n2 component=time-vs-fitness dataset=h3n2_clades static=figures/h3n2_clades_time_vs_fitness.png}
**Cumulative H3N2 fitness flux.**
Empirical frequencies of H3N2 clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

SARS-CoV-2 clade frequencies and fitnesses can be viewed as a phase portrait, plotting each clade's empirical frequency against its fitness relative to the daily population average ([@fig:frequency-vs-fitness-sarscov2]).
A clade emerges at low frequency and high relative fitness, sweeps up in frequency as its relative fitness declines toward the population average, peaks near a relative fitness of zero, and then falls back to low frequency as it is outcompeted.
Clades that start out with a greater advantage over the population average tend to sweep to higher maximum frequency than clades that start with less of an advantage.

:::figure{#fig:frequency-vs-fitness-sarscov2 component=frequency-vs-fitness dataset=sarscov2_clades}
**Frequency vs fitness phase diagram for SARS-CoV-2 clades.**
Each line traces a SARS-CoV-2 clade's trajectory over time through empirical frequency (x-axis, logit scale) and fitness relative to the daily population average (y-axis), estimated from multinomial logistic regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

### Fisher's fundamental theorem

Multistrain models that allow for antigenic evolution produce traveling waves in antigenic space [@bedford2012canalization].
More broadly, many mutations of small fitness effect create traveling fitness waves where the the rate of translation to higher fitness is proportional to the variance in fitness [@neher2013genealogies].
This is a consequence of Fisher's fundamental theorem of natural selection
$$\frac{d\bar{f}}{dt} = \mathrm{Var}(f),$$
where "the rate of increase in fitness of any organism at any time is equal to its genetic variance in fitness at that time" [@fisher1930genetical].

We can investigate this relationship directly in SARS-CoV-2 ([@fig:sarscov2-variance-flux]), where we find that timepoints with larger variance in fitness $\mathrm{Var}[f(t)] = \sum_i x_i(t) \, (f_i - \bar{f}(t))^2$ correlate well with timepoints with larger change in mean population fitness $\Delta \bar{f}(t) / \Delta t$.
In fact we find that the relationship is near the 1:1 expectation from Fisher's theorem  (slope = 1.20, Pearson $r$ = 0.96).
Looking in detail at rate of fitness flux through time, we find a peak fitness flux of $8.5 \times 10^{-3}$ per-gen in 2021 followed by a reduction to $1.6-1.7 \times 10^{-3}$ per-gen in 2024 and 2025.
This shows that the rate of adaptation of SARS-CoV-2 has been slowing as low hanging fruit of host adaptation is exhausted, leaving only red-queen antigenic evolution to drive adaptation.

:::figure{#fig:sarscov2-variance-flux component=variance-vs-flux dataset=sarscov2_clades scalemax=40 static=figures/sarscov2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in SARS-CoV-2.**
Fitness variance is compared to fitness flux, where each dot represents a daily timepoint.
:::

Compared to SARS-CoV-2, influenza H3N2 shows generally lower rates of fitness flux, averaging $0.5 \times 10^{-3}$ per-gen from 2016 to 2026 ([@fig:h3n2-variance-flux]).
This is roughly 3 times lower than recent years of SARS-CoV-2 fitness flux.
However, it remains possible that SARS-CoV-2 slows further in the coming years.

:::figure{#fig:h3n2-variance-flux component=variance-vs-flux dataset=h3n2_clades static=figures/h3n2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in H3N2.**
Fitness variance is compared to fitness flux, where each dot represents a daily timepoint.
:::

For both SARS-CoV-2 and influenza H3N2, the connection between fitness flux and fitness variance is clear.
This suggests from first principles that interventions that decrease variance in fitness across the virus population would be expected to slow adaptation, while interventions that increase variance would be expected to speed adaptation.

### Mutational fitness effects

A simple correlation of cumulative mutations against cumulative fitness flux will fail due to phylogenetic non-independence, and instead we rely on phylogenetic contrasts of parent and daughter lineages [@felsenstein1985phylogenies].
Pango lineages [@rambaut2020dynamic] provide a convenient granular and hierarchical nomenclature well suited to this.

We define Pango parent-to-child branches by finding mutations that accumulate between hierarchical Pango lineages.
Given estimates of per-lineage fitness, we also calculate the change in fitness between parent/child pairs.
As expected, we observe that Pango lineages are granular with parent-to-child changes in spike mutations, non-spike mutations and fitness being modest, with most branches adding only a handful of substitutions and shifting fitness by a small amount ([@fig:delta-hist]).

:::figure{#fig:delta-hist component=lineage-delta-histograms dataset=sarscov2_lineages predictors=spike,nonspike}
**Distributions of mutation and fitness change across SARS-CoV-2 lineage branches.**
Across all parent-to-child Pango lineage branches, the change in the number of spike substitutions, the change in non-spike substitutions, and the change in fitness.
Each bar is the fraction of branches in that bin; rare large founder jumps fall beyond the plotted range.
:::

We compare the change in fitness along each parent-to-child branch against the substitutions it acquired in different regions of the SARS-CoV-2 genome ([@fig:delta-genome]).
Spike S1 substitutions and those in the receptor-binding domain (RBD) in particular, carry the strongest positive association with fitness gain, while ORF1ab and accessory gene substitutions show lower regression slopes and correlation coefficients.
Splitting branches into early and late periods shows that the fitness value of a given substitution is largest early and decays over time.

:::figure{#fig:delta-genome component=lineage-deltas dataset=sarscov2_lineages predictors=s1,rbd,orf1ab,accessory}
**Lineage-specific amino acid change versus lineage-specific fitness change across regions of the SARS-CoV-2 genome.**
Each point is one parent-to-child Pango lineage branch in one season: the change in the number of substitutions in a genome region (x) against the change in fitness (y), colored by time from blue (2020) to red (2026), with a least-squares fit per panel.
The All / Early / Late toggle restricts to early (Jan 2020–Jun 2022) or late (Jul 2022 onward) branches.
:::

We can track this relationship through time by fitting the regression separately within each season ([@fig:delta-trends]).
The per-substitution effect on fitness of spike and RBD changes is largest early and erodes toward zero as the readily accessible routes to host adaptation are exhausted, while other regions stay near zero throughout.
However, note that across accessory proteins and in ORF1ab there is moderate correlation from 2020 to 2022 of mutations to changes in fitness suggesting some lesser, but likely non-zero, involvement in adaptation compared to spike.
Although per-substitution effect (ie regression slope) of RBD decays from 2020, the predictive ability of spike S1 and RBD substitutions as measured by Pearson and Spearman correlations stays high through the period with average correlation coefficients of $r$ = 0.73 and $\rho$ = 0.56 for spike S1 and $r$ = 0.68 and $\rho$ = 0.51 for spike RBD.

:::figure{#fig:delta-trends component=lineage-delta-trends dataset=sarscov2_lineages predictors=s1,rbd,orf1ab,accessory}
**Strength of the mutation–fitness relationship through time.**
For each season the parent-to-child branches are summarized into one statistic relating change in regional substitutions to change in fitness, with one line per genome region.
Toggle between the regression slope, Pearson *r*, and Spearman *ρ*.
:::

The marginal relationships in [@fig:delta-genome] cannot on their own establish which substitutions drive fitness.
A more evolved lineage more accumulates substitutions across the whole genome, so a region can correlate with fitness merely by co-varying with a genuinely causal region.
The moderate marginal association of ORF1ab substitutions is a case in point.
To isolate each region's independent contribution we fit a multiple linear regression of the change in substitution count across four non-overlapping genome regions to per-branch change in fitness ([@fig:delta-multiple]).
Once spike is controlled for, essentially all of the positive signal sits in the RBD ($\beta$ = 0.055 per substitution) with a smaller contribution from the remainder of S1 ($\beta$ = 0.016), while the ORF1ab coefficient collapses to near zero and is no longer distinguishable from no effect ($\beta$ = 0.002, $p$ = 0.24), as is the accessory-protein coefficient ($\beta$ = –0.001, $p$ = 0.48).
The apparent marginal association of ORF1ab is therefore a confound of its co-occurrence with spike change rather than evidence that ORF1ab substitutions themselves raise fitness.
The four-region model explains a majority of the variance in per-branch fitness change ($R^2$ = 0.59), with predicted and observed changes falling along the 1:1 line.

:::figure{#fig:delta-multiple component=lineage-deltas-model dataset=sarscov2_lineages}
**Multiple regression of fitness change on non-overlapping genome regions.**
A pooled ordinary-least-squares fit of the per-branch change in fitness on the change in substitution count in four non-overlapping regions of the SARS-CoV-2 genome (spike RBD, spike S1 outside the RBD, ORF1ab, accessory).
The table gives each region's partial estimate.
The scatter plots each parent-to-child Pango lineage branch's model-predicted change in fitness (x) against its observed change in fitness (y), colored by time from blue (2020) to red (2026) and with dashed 1:1 calibration line.
The All / Early / Late toggle refits the model over all branches or the early (Jan 2020–Jun 2022) and late (Jul 2022 onward) subsets, updating both the table and the scatter.
:::

### Predicting fitness effects

The core idea of comparing change in mutation count to change in fitness expresses a similar logic to McDonald-Kreitman tests [@mcdonald1991adaptive] where the key comparison is relative success of lineages bearing different mutation patterns.
This predictor vs growth rate formulation [@kistler2022rapid] should be robust to many sources of confounding that other metrics testing drivers of adaptive evolution suffer from.

We come our simple non-parameterized fitness predictor of relative spike S1 substitution count to recent proposals to predict evolutionary successful substitutions via deep learning ([@fig:delta-predictors]).
EvEscape [@thadani2023learning] uses a variational autoencoder alongside accessibility and biochemical dissimilarity to score SARS-CoV-2 spike mutations. 
Semanticity measures Euclidean distance of protein language model embeddings [@hie2021learning].
Here, we use pre-computed per-lineage EvEscape score alongside a reimplemention of semanticity using ESM-2 protein language model embeddings [@lin2023evolutionary].
We analyze both the pre-trained 650M parameter ESM-2 model as well as model fine-tuned to 16k SARS-CoV-2 spike sequences with collection dates from 2020 through 2022.
We find that both EvEscape score as well as fine-tuned ESM-2 semanticity perform similarly to simple count of spike S1 substitutions to disambiguate fitness of circulating SARS-CoV-2 lineages.

:::figure{#fig:delta-predictors component=lineage-deltas dataset=sarscov2_lineages predictors=s1,evescape,esm_650M_pretrained,esm_650M_fine_tuned}
**Lineage-specific predictors versus lineage-specific fitness change.**
Each point is one parent-to-child Pango lineage branch in one season: change in predictor value (x) against the change in fitness (y), colored by time from blue (2020) to red (2025), with a least-squares fit per panel.
The All / Early / Late toggle restricts to early (Jan 2020–Jun 2022) or late (Jul 2022 onward) branches.
:::

## Conclusions

Most measures of viral adaptation are indirect, diagnosing the presence of selection based on mutations patterns.
Here we instead read adaptation directly off the dynamics of variant frequencies, aggregating per-variant growth rates into the population's fitness flux.
This turns the tempo of adaptation into a single quantity that can be followed through time, placed on a common per-generation scale across pathogens and connected to first principles.
On this scale SARS-CoV-2 adapts rapidly, doubling in fitness roughly every 1.5 years, but decelerating from a peak flux in 2021 towards a baseline flux in 2024, while seasonal H3N2 sustains a slower, steadier flux.
Importantly these numbers can be connected back to epidemiological impacts [@figgins2025frequency] and have an absolute scale to them.

The generality of the approach rests on a single requirement: a way to bin genetic diversity into discrete, comparable variants.
For SARS-CoV-2 and influenza this comes off the shelf, with Nextstrain clades supporting the frequency and flux analysis and finer, hierarchically nested Pango lineages supporting the phylogenetic-contrast analysis of mutational effects.
Pathogens without an established nomenclature could be analyzed via automated methods that partition a tree into lineages [@mcbroome2024framework; @lefrancq2025learning].

Beyond describing historical adaptation, the per-branch contrast of mutation against fitness change yields a simple, interpretable account of which substitutions matter.
For SARS-CoV-2 the signal concentrates in spike (and particularly spike RBD) and a plain count of spike S1 substitutions disambiguates the relative fitness of co-circulating lineages about as well as recently proposed deep-learning escape and protein-language-model scores [@thadani2023learning; @hie2021learning].
This makes the mutation-to-fitness deltas a strong and transparent baseline for forecasting variant success, where a predictor that does not improve on counting spike substitutions has not yet justified its added complexity.

## Methods

### Sequence data

For SARS-CoV-2, we use curated "open" data from Nextstrain [@hadfield2018nextstrain] that draws from NCBI GenBank.
For influenza H3N2, we use data from GISAID [@shu2017gisaid].
In each case, the raw sequences are processed with Nextclade [@aksamentov2021nextclade] to assign Nextstrain clade and Pango lineage [@rambaut2020dynamic] to SARS-CoV-2 sequences and to assign subclade [@neher2026nomenclature] to influenza H3N2 sequences.
We filter out sequences with Nextclade overall QC status of "bad".
We additionally filter to sequences collected from the USA.
This leaves 3,594,555 total sequences for SARS-CoV-2 sampled between 2020 and 2026 and 49,623 total sequences for H3N2 sampled between 2016 and 2026.

We conducted multinomial logistic regression (MLR) using the evofr package ([github.com/blab/evofr](https://github.com/blab/evofr)) on 1-year sliding windows for SARS-CoV-2 (12 windows total) and 2-year sliding windows for H3N2 (10 windows total).
For each window we treat each clade as a distinct variant, collapsing rare clades together into a single "other" category before fitting.
For both SARS-CoV-2 and H3N2, a clade is modeled separately only if it reaches at least 50 sequences and a mean frequency of at least 0.1% within the window, while clades below either threshold are merged into "other".
This leaves between 7 and 18 clades per window (median 15) for SARS-CoV-2 and between 5 and 13 (median 9) for influenza H3N2.

We conduct a parallel MLR analysis of SARS-CoV-2 Pango lineages.
Because lineages are hierarchically nested, rather than collapsing rare lineages into a shared "other" we roll each lineage with fewer than 500 sequences up into its parent lineage, repeating until every retained lineage clears this count.
A lineage is additionally retained only if at least 200 sequences are assigned to that lineage itself rather than to a descendant sub-lineage, otherwise it is folded into "other".
This leaves between 13 and 165 lineages per window (median 79) for SARS-CoV-2.
Rationale for specific collapse cutoffs is available at [github.com/blab/fitness-flux/tree/main/inclusion-thresholds](https://github.com/blab/fitness-flux/tree/main/inclusion-thresholds).

### Scaffolding across timepoints

Within each sliding window the MLR model estimates each variant's fitness only relative to that window's pivot, so every window carries its own arbitrary additive zero and the per-window estimates $f_{i,w}$ are not directly comparable.
We recover a single fitness per variant by treating scaffolding as a weighted two-way additive model: each estimate is a variant effect minus a window effect, $f_{i,w} \approx f_i - c_w$, where $f_i$ is variant $i$'s global fitness and $c_w$ is window $w$'s offset.
We choose the $f_i$ and $c_w$ that jointly minimize the abundance-weighted squared error across every window,
$$\min_{\{f_i\},\,\{c_w\}} \; \sum_{i,w} a_{i,w} \, (f_{i,w} - f_i + c_w)^2,$$
weighting each estimate by $a_{i,w}$, the area under variant $i$'s modeled-frequency curve in window $w$, so a window in which a variant is rare with poorly constrained MLR estimate will contribute negligibly.
The optimum is a pair of interleaved abundance-weighted means,
$$f_i = \frac{\sum_w a_{i,w} \, (f_{i,w} + c_w)}{\sum_w a_{i,w}}, \qquad c_w = \frac{\sum_i a_{i,w} \, (f_i - f_{i,w})}{\sum_i a_{i,w}},$$
each variant's fitness being the weighted mean of its offset-corrected estimates over its windows, and each window's offset the weighted-mean gap between the global scale and that window's estimates; we solve by alternating the two to convergence.
The overlap of variants between windows ties them into one connected scale, leaving a single global constant free, which we fix by shifting all values so that the founding variant, our least-fit baseline, sits at zero, leaving each variant's scaffolded value as its cumulative fitness flux $\Phi_i = f_i - f_0$.

### Lineage mutation counts and branch contrasts

To relate change in fitness to change in genotype, we count amino-acid substitutions per Pango lineage and compare each lineage against its parent.
Per-lineage substitution counts are read from the Nextclade SARS-CoV-2 reference tree at [nextstrain.org/nextclade/nextstrain/sars-cov-2/wuhan-hu-1/orfs](https://nextstrain.org/nextclade/nextstrain/sars-cov-2/wuhan-hu-1/orfs), in which each tip corresponds to a Pango lineage.
For each lineage we count substitutions relative to the Wuhan-Hu-1 reference.
We tally these by region separating: spike S1 subunit, the receptor-binding domain (RBD, 319–541) within S1, ORF1ab and the accessory and structural genes (ORF3a, E, M, ORF6, ORF7a, ORF7b, ORF8, N).

We then form parent-to-child branches between hierarchically nested Pango lineages.
Within each window a lineage's parent is its closest retained ancestor, where the retained set is fixed by the collapsing described above, so lineages that were rolled up into a parent or folded into "other" do not themselves appear as branch endpoints.
For every branch whose parent and child both carry an MLR fitness estimate in that window, we record the change in substitution count in each genome region and the change in fitness, taken as the difference in their per-window MLR fitness; because both endpoints are estimated against the same window's pivot, this contrast is well defined without scaffolding.
A branch is recorded once per window in which it appears, so a lineage pair that co-circulates across several windows contributes several observations.
These per-branch deltas in mutation count and fitness are the unit of the mutational-fitness analyses.

### Mutational fitness predictors

We compare the per-branch substitution counts above against two externally proposed predictors of mutational fitness, each reduced to a per-lineage value and contrasted across the same parent-to-child branches. EvEscape [@thadani2023learning] combines a variational-autoencoder fitness model with residue accessibility and biochemical dissimilarity to score spike mutations.
We use the precomputed all-strain EvEscape scores released at [evescape.org/data](https://evescape.org/data) and take each Pango lineage's score as the mean EvEscape across the sequences assigned to that lineage.

Semanticity follows the semantic-change measure of Hie et al. [@hie2021learning], reimplemented with ESM-2 [@lin2023evolutionary].
Each Pango lineage's spike amino-acid sequence is embedded with the 650M-parameter ESM-2 model (`esm2_t33_650M_UR50D`), taking the CLS-token representation of the final (33rd) layer as a 1280-dimensional sequence embedding.
We embed lineages with both the released pretrained weights and weights fine-tuned under a masked-language-model objective (15% of residues masked) for one epoch (AdamW, learning rate $5 \times 10^{-5}$) on roughly 16,000 SARS-CoV-2 spike sequences collected from 2020 through 2022.
Fine-tuning and embedding code is available at [github.com/blab/embedded-pathways](https://github.com/blab/embedded-pathways).
The semanticity of a branch is the Euclidean distance between its child and parent lineage embeddings.

## Acknowledgments

SARS-CoV-2 analyses are based on open data in GenBank.
We gratefully acknowledge the researchers and data contributors who collected the specimens, generated and deposited the raw sequence data and metadata into NCBI GenBank.
Influenza analyses are based on GISAID data.
We gratefully acknowledge all data contributors, i.e., the Authors and their Originating laboratories
responsible for obtaining the specimens, and their Submitting laboratories for generating the
genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is
based. 
TB was funded as a Howard Hughes Medical Institute Investigator.
This research was supported in part by grant NSF PHY-2309135, the Gordon and Betty Moore Foundation grant no. 2919.02, and the Chan Zuckerberg Initiative DAF grant to the Kavli Institute for Theoretical Physics (KITP).
