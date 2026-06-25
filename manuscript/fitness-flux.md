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

TBD

## Introduction

RNA viruses evolve rapidly, and the selective pressures they face shift over the course of emergence to endemicity.
A virus newly established in humans initially adapts to its new host by refining its capacity for replication and transmission.
Once a virus becomes endemic, adaptation in some viruses is instead dominated by continual escape from accumulating population immunity, driving ongoing antigenic change.
Adaptation of either kind leaves a signature in nonsynonymous vs synonymous substitutions, with methods ranging from simple comparisons of nonsynonymous to synonymous substitution rates (dN/dS) to McDonald–Kreitman-style approaches that weigh mutations fixed along a virus's successful trunk lineage against those lost on unsuccessful side branches [@wolf2006long].
Such approaches have revealed rapid, parallel adaptation in the SARS-CoV-2 spike S1 subunit [@kistler2022rapid] and continuous adaptive evolution across endemic human viruses more broadly [@kistler2023atlas].

A complementary class of methods estimates fitness directly from the dynamics of variant frequencies rather than from the composition of mutations.
Multinomial logistic regression (MLR) models the frequencies of co-circulating variants through time and infers a relative growth rate, or fitness, for each [@obermeyer2022analysis; @abousamra2024fitness].
Because it expresses fitness as a difference in growth rate between variants, this measure maps directly onto the population-genetic notion of selective advantage.
These growth-rate differences correspond to differences in the time-varying effective reproduction number between co-circulating variants [@figgins2025frequency].
Aggregating these per-variant fitnesses into the rate of change of mean population fitness yields the population's fitness flux [@mustonen2010fitness], a direct, frequency-based alternative to dN/dS for quantifying the tempo of adaptation.

Here we use this frequency-based view of fitness to trace how SARS-CoV-2 has adapted from the early pandemic in 2020 through to the present, spanning the transition from initial host adaptation to sustained evolution for antigenic novelty.
We place the rate of SARS-CoV-2 fitness change in context by comparing it against seasonal influenza A/H3N2, a canonical rapidly and continuously adapting human virus.
Finally, we relate the inferred changes in fitness to molecular predictors, most directly the accumulation of spike mutations, to identify the substitutions that drive fitness gain.

## Results and discussion

### Frequency dynamics

We follow population genetics first principles to compute the frequency through time of a haploid allele under selection.
If an allele is at frequency $x$, then after a single generation with selective advantage $s$ the expected allele frequency will be 
$$x' = \frac{x \, (1+s)}{x \, (1+s) + (1-x)}.$$
Compounded over $t$ generations, the expectation from initial frequency $p$ follows
$$x(t) = \frac{p \, (1+s)^t}{p \, (1+s)^t + (1-p)}.$$
Generalizing this two-type model to $n$ co-circulating variants, each with initial frequency $p_i$ and selective advantage $s_i$, variant $i$'s frequency is its relative size normalized by the sum across all variants,
$$x_i(t) = \frac{p_i \, (1+s_i)^t}{\sum_j p_j \, (1+s_j)^t},$$
where the two-type model above is recovered by a single focal allele competing against a reference type with $s=0$.
Moving from discrete generations to continuous time, $(1+s_i)^t = \mathrm{exp}(t \, \mathrm{log}(1+s_i))$, so writing the growth rate $f_i = \mathrm{log}(1+s_i)$ gives the probability that a virus sampled at time $t$ is labeled as variant $i$
$$\mathrm{Pr}(X = i) = x_i(t) = \frac{p_i \, \mathrm{exp}(f_i \, t)}{\sum_j p_j \, \mathrm{exp}(f_j \, t) }.$$

This is the multinomial logistic regression (MLR) model, which has seen significant previous use for modeling SARS-CoV-2 variant frequencies [@abousamra2024fitness].
The denominator normalizes the exponential growth/decay of individual variants so that overall frequency sums to 1, and the model has $2n$ parameters, with each variant $i$ having an initial frequency $p_i$ and a fixed growth rate $f_i$.
Because growth rates are necessarily relative, we fix an arbitrary "pivot" variant as a reference with growth rate $f=0$.

We estimate frequencies and fitnesses of SARS-CoV-2 clades in 1-year sliding windows between Jan 2020 and Jan 2026 ([@fig:time-vs-frequency-sarscov2]).
In each window we collect sequence counts for SARS-CoV-2 clades from sequences from the USA and estimate per-variant frequencies and fitnesses.
The match between the empirical frequencies (dotted trajectories) and MLR frequencies (solid trajectories) indicates the model fits well despite having few parameters.

:::figure{#fig:time-vs-frequency-sarscov2 component=time-vs-frequency dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_frequency.png}
**Relative frequencies of SARS-CoV-2 clades through time.**
Points represent empirical frequencies of SARS-CoV-2 Nextstrain clades, while solid lines represent modeled frequencies from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

We take a similar approach to estimating frequencies and fitnesses of seasonal influenza H3N2 ([@fig:time-vs-frequency-h3n2]).
However, here we use 2-year sliding windows to account for slower frequency dynamics in seasonal influenza.
The model fits are not as good for H3N2 compared to SARS-CoV-2.
This is especially apparent at junctions between influenza seasons where stochastic seeding of new season may result in a discontinuity of clade frequency compared to MLR expectation.

:::figure{#fig:time-vs-frequency-h3n2 component=time-vs-frequency dataset=h3n2_clades static=figures/h3n2_clades_time_vs_frequency.png}
**Relative frequencies of H3N2 clades through time.**
Points represent empirical frequencies of H3N2 Nextstrain clades, while solid lines represent modeled frequencies from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

We "scaffold" MLR fitness estimates across windows to arrive at a single fitness estimate per variant.
Here, each window only measures fitness differences, with its own arbitrary zero; we solve for the single set of clade fitnesses and per-window offsets that best fit every window at once, weighting by abundance (see Methods).

### Fitness flux

With variant frequency $x_i(t)$ and constant variant fitness $f_i$, we can describe the mean population fitness as a standard weighted sum $\bar{f}(t) = \sum_i x_i(t) \, f_i$.
The fitness flux [@mustonen2010fitness] of the population is then the rate of change of population fitness at a given time $\phi(t) = \Delta \bar{f}(t) / \Delta t$.
Integrating this rate gives the cumulative fitness flux
$$\Phi(t) = \int_{t_0}^{t} \phi(t') \, dt' = \bar{f}(t) - \bar{f}(t_0),$$
the total adaptive change accumulated along the population's trajectory.
Because variant fitnesses are estimated only relative to a pivot, an individual variant's scaffolded fitness is meaningful as a difference from a baseline rather than as an absolute value.
Chaining these locally-measured advantages across overlapping windows places variant $i$ at a cumulative fitness flux $\Phi_i = f_i - f_0$ relative to the founding variant, and the population sits at the frequency-weighted average $\Phi(t) = \sum_i x_i(t) \, \Phi_i$.

We see SARS-CoV-2 accumulates fitness flux rapidly, over the course of 2020 through 2025, doubling in fitness every 1.5 years ([@fig:time-vs-fitness-sarscov2]).
We see an initial lull, followed by rapid growth in fitness in 2021 and 2022, and then a slower, more steady pace since 2024.
There is a mix of large jumps in fitness (in particular to BA.1, but also more recently to JN.1) and smaller, more gradual step change.

:::figure{#fig:time-vs-fitness-sarscov2 component=time-vs-fitness dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_fitness.png}
**Cumulative SARS-CoV-2 fitness flux.**
Empirical frequencies of SARS-CoV-2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

Seasonal influenza H3N2 shows a fundamentally similar pattern of emergence of new clades and replacement of existing diversity, however, it plays out on a slower timescale ([@fig:time-vs-fitness-h3n2]).
Rather than doubling in fitness every 1.5 years, H3N2 has a doubling period of 9.6 years.
More coexistence of multiple co-circulating clades is also apparent relative to SARS-CoV-2.

:::figure{#fig:time-vs-fitness-h3n2 component=time-vs-fitness dataset=h3n2_clades static=figures/h3n2_clades_time_vs_fitness.png}
**Cumulative H3N2 fitness flux.**
Empirical frequencies of H3N2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

SARS-CoV-2 frequencies and fitnesses can be viewed as a phase portrait, plotting each clade's empirical frequency against its fitness relative to the daily population average ([@fig:frequency-vs-fitness-sarscov2]).
A clade emerges at low frequency and high relative fitness, sweeps up in frequency as its relative fitness declines toward the population average, peaks near a relative fitness of zero, and then falls back to low frequency as it is outcompeted.
Clades that start out with a greater advantage over the population average tend to sweep to higher maximum frequency than clades that start with less of an advantage.

:::figure{#fig:frequency-vs-fitness-sarscov2 component=frequency-vs-fitness dataset=sarscov2_clades}
**Frequency vs fitness phase diagram for SARS-CoV-2 clades.**
Each line traces a SARS-CoV-2 Nextstrain clade's trajectory over time through empirical frequency (x-axis, logit scale) and fitness relative to the daily population average (y-axis), estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

### Fisher's fundamental theorem

Multistrain models that allow for antigenic evolution produce traveling waves in antigenic space [@bedford2012canalization].
More broadly, many mutations of small fitness effect create traveling fitness waves where the the rate of translation to higher fitness is proportional to the variance in fitness [@neher2013genealogies].
This is a consequence of Fisher's fundamental theorem of natural selection
$$\frac{d\bar{f}}{dt} = \mathrm{Var}(f),$$
or "the rate of increase in fitness of any organism at any time is equal to its genetic variance in fitness at that time" [@fisher1930genetical].

If we investigate this relationship directly in SARS-CoV-2 ([@fig:sarscov2-variance-flux]), we find that timepoints with larger variance in fitness $\mathrm{Var}[f(t)] = \sum_i x_i(t) \, (f_i - \bar{f}(t))^2$ correlate well with timepoints with larger change in mean population fitness $\Delta \bar{f}(t) / \Delta t$.
In fact we find that the relationship is near the 1:1 expectation from Fisher's theorem.
Looking in detail at rate of fitness flux through time, we see reduction to per gen average of $1.6-1.7 \times 10^{-3}$ in 2024 and 2025, down from $8.5 \times 10^{-3}$ in 2021.
This shows that the rate of adaptation of SARS-CoV-2 has been slowing as low hanging fruit of host adaptation is exhausted, leaving only red-queen antigenic evolution to drive adaptation.

:::figure{#fig:sarscov2-variance-flux component=variance-vs-flux dataset=sarscov2_clades scalemax=40 static=figures/sarscov2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in SARS-CoV-2.**
Log fitness variance is compared to log fitness flux, where each dot represents a daily timepoint.
:::

Compared to SARS-CoV-2, influenza H3N2 shows generally lower rates of fitness flux, averaging $0.5 \times 10^{-3}$ from 2016 to 2026 ([@fig:h3n2-variance-flux]).
This is roughly 3 times lower than recent years of SARS-CoV-2 fitness flux, but it's possible that SARS-CoV-2 slows further.

:::figure{#fig:h3n2-variance-flux component=variance-vs-flux dataset=h3n2_clades static=figures/h3n2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in H3N2.**
Log fitness variance is compared to log fitness flux, where each dot represents a daily timepoint.
:::

For both SARS-CoV-2 and influenza H3N2, the connection between fitness flux and fitness variance is clear.
This suggests from first principles that interventions that decrease variance in fitness across the virus population would be expected to slow adaptation, while interventions that increase variance would be expected to speed adaptation.

### Mutational fitness effects

We shouldn't do a simple correlation of mutations against cumulative fitness flux due to phylogenetic non-independence.
Instead can rely on phylogenetic contrasts of parent and daughter lineages.
Pango lineages provide a convenient granular and hierarchical nomenclature well suited to this.

Each parent-to-child branch contributes a change in mutation count and a change in fitness; across all branches these are modest, with most branches adding only a handful of substitutions and shifting fitness by a small amount ([@fig:delta-hist]).

:::figure{#fig:delta-hist component=lineage-delta-histograms dataset=sarscov2_lineages predictors=spike,nonspike}
**Distributions of mutation and fitness change across SARS-CoV-2 lineage branches.**
Across all parent-to-child Pango lineage branches, the change in the number of spike substitutions, the change in non-spike substitutions, and the change in log fitness.
Each bar is the fraction of branches in that bin; rare large founder jumps fall beyond the plotted range.
:::

We compare the change in fitness along each parent-to-child lineage branch against the substitutions it acquired in different regions of the SARS-CoV-2 genome ([@fig:delta-genome]).
Spike substitutions, and those in the receptor-binding domain (RBD) in particular, carry the strongest positive association with fitness gain, while accessory-gene substitutions carry almost none.
Splitting branches into early and late periods shows that the fitness value of a given substitution is largest early and erodes over time as the population approaches the limits of antigenic escape.

:::figure{#fig:delta-genome component=lineage-deltas dataset=sarscov2_lineages predictors=s1,rbd,orf1ab,accessory}
**Lineage-specific amino acid change versus lineage-specific fitness change across regions of the SARS-CoV-2 genome.**
Each point is one parent-to-child Pango lineage branch in one season: the change in the number of substitutions in a genome region (x) against the change in log fitness (y), colored by time from blue (2020) to red (2025), with a least-squares fit per panel.
The All / Early / Late toggle restricts to early (Jan 2020–Jun 2022) or late (Jul 2022 onward) branches.
:::

We can track this relationship through time by fitting the regression separately within each season ([@fig:delta-trends]).
The fitness value of spike and RBD substitutions is largest early and erodes toward zero as the readily accessible routes to host adaptation are exhausted, while other regions stay near zero throughout.
However, note that across accessory proteins and in ORF1ab there is moderate correlation from 2020 to 2022 of mutations to changes in fitness suggesting some lesser, but likely non-zero, involvement in adaptation compared to spike.

:::figure{#fig:delta-trends component=lineage-delta-trends dataset=sarscov2_lineages predictors=s1,rbd,orf1ab,accessory}
**Strength of the mutation–fitness relationship through time.**
For each season the parent-to-child branches are summarized into one statistic relating change in regional substitutions to change in log fitness, with one line per genome region.
Toggle between the regression slope, Pearson *r*, and Spearman *ρ*; computed over the same branches as [@fig:delta-genome].
:::

### Predicting fitness effects

:::figure{#fig:delta-predictors component=lineage-deltas dataset=sarscov2_lineages predictors=s1,evescape,esm_650M_pretrained,esm_650M_fine_tuned}
**Lineage-specific predictors versus lineage-specific fitness change.**
Each point is one parent-to-child Pango lineage branch in one season: change in predictor value (x) against the change in log fitness (y), colored by time from blue (2020) to red (2025), with a least-squares fit per panel.
The All / Early / Late toggle restricts to early (Jan 2020–Jun 2022) or late (Jul 2022 onward) branches.
:::

## Conclusions

TBD

## Methods

### Sequence data

For SARS-CoV-2 we use curated "open" data from Nextstrain [@hadfield2018nextstrain] that draws from NCBI GenBank.

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

## Acknowledgments

This research was supported in part by grant NSF PHY-2309135, the Gordon and Betty Moore Foundation grant no. 2919.02, and the Chan Zuckerberg Initiative DAF grant to the Kavli Institute for Theoretical Physics (KITP).
