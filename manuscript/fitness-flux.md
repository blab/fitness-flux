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

TBD

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

With variant frequency $x_i(t)$ and constant variant fitness $f_i$, we can describe the mean population fitness as a standard weighted sum $\bar{f}(t) = \sum_i x_i(t) \, f_i$.
The fitness flux [@mustonen2010fitness] of the population is then the rate of change of population fitness at a given time $\phi(t) = \Delta \bar{f}(t) / \Delta t$.



:::figure{#fig:time-vs-fitness-sarscov2 component=time-vs-fitness dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_fitness.png}
**Cumulative SARS-CoV-2 fitness flux.**
Empirical frequencies of SARS-CoV-2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::


:::figure{#fig:time-vs-fitness-h3n2 component=time-vs-fitness dataset=h3n2_clades static=figures/h3n2_clades_time_vs_fitness.png}
**Cumulative H3N2 fitness flux.**
Empirical frequencies of H3N2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

The same frequencies and fitnesses can be viewed as a phase portrait, plotting each clade's empirical frequency against its fitness relative to the daily population average ([@fig:frequency-vs-fitness-sarscov2]).
A clade emerges at low frequency and high relative fitness, sweeps up in frequency as its relative fitness declines toward the population average, peaks near a relative fitness of zero, and then falls back to low frequency as it is outcompeted.
Clades that start out with a greater advantage over the population average tend to sweep to higher maximum frequency than clades that start with less of an advantage.

:::figure{#fig:frequency-vs-fitness-sarscov2 component=frequency-vs-fitness dataset=sarscov2_clades}
**Frequency vs fitness phase diagram for SARS-CoV-2 clades.**
Each line traces a SARS-CoV-2 Nextstrain clade's trajectory over time through empirical frequency (x-axis, logit scale) and fitness relative to the daily population average (y-axis), estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

Fisher's fundamental theorem and the expected relationship of fitness variance and fitness flux.
Empirical investigation of this relationship in SARS-CoV-2 and H3N2.

:::figure{#fig:sarscov2-variance-flux component=variance-vs-flux dataset=sarscov2_clades static=figures/sarscov2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in SARS-CoV-2.**
Log fitness variance is compared to log fitness flux, where each dot represents a daily timepoint.
:::

:::figure{#fig:h3n2-variance-flux component=variance-vs-flux dataset=h3n2_clades static=figures/h3n2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in H3N2.**
Log fitness variance is compared to log fitness flux, where each dot represents a daily timepoint.
:::

We shouldn't do a simple correlation of mutations against cumulative fitness flux due to phylogenetic non-independence.
Instead can rely on phylogenetic contrasts of parent and daughter lineages.
Pango lineages provide a convenient granular and hierarchical nomenclature well suited to this.

:::figure{#fig:delta-genome src=figures/sarscov2_lineage_delta_fitness_across_genome.png}
**Correlation of lineage-specific amino acid change to lineage-specific fitness change across regions of the SARS-CoV-2 genome.**
:::

:::figure{#fig:delta-time src=figures/sarscov2_lineage_delta_fitness_across_time.png}
**Correlation of lineage-specific amino acid change to lineage-specific fitness change over time.**
:::

:::figure{#fig:delta-evescape src=figures/sarscov2_lineage_delta_fitness_vs_evescape.png}
**Correlation of lineage-specific change in EvEscape score to lineage-specific fitness change.**
:::

## Conclusions

TBD

## Methods

TBD

## Acknowledgments

This research was supported in part by grant NSF PHY-2309135, the Gordon and Betty Moore Foundation grant no. 2919.02, and the Chan Zuckerberg Initiative DAF grant to the Kavli Institute for Theoretical Physics (KITP).
