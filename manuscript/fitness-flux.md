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
If the allele is initial at frequency $x$ after a single generation with selective advantage $s$ the expected allele frequency will be 
$$x' = \frac{x \, (1+s)}{x \, (1+s) + (1-x)}.$$
Compounded over $t$ generations, the expectation from initial frequency $x_0$ follows
$$x(t) = \frac{x_0 \, (1+s)^t}{x_0 \, (1+s)^t + (1-x_0)}.$$
Here, we can see that trajectories are linear once logit transformed via $\mathrm{log}(\frac{x}{1 - x})$.

:::figure{#fig-time-vs-frequency-sarscov2 src=figures/pop_gen_logit_trajectories.png}
**Expected trajectories of a selected allele in a haploid population genetics model.**
The left hand panel shows with variant frequency $x$ in normal-space, while the right-hand panel shows logit-transformed variant frequency.
Selective advantage $s$ of 1%, 2% and 5% per generation are shown.
:::

We follow this population genetics logic in implementing multinomial logistic regression (MLR), which has seen significant previous use for modeling SARS-CoV-2 variant frequencies [@abousamra2024fitness].
MLR across $n$ variants models the probability of a virus sampled at time $t$ to be labeled as variant $i$ as equal to its frequency $x_i(t)$
$$\mathrm{Pr}(X = i) = x_i(t) = \frac{p_i \, \mathrm{exp}(f_i \, t)}{\sum_j p_j \, \mathrm{exp}(f_j \, t) },$$
where the denominator serves to normalize exponential growth/decay of individual variants and keep overall frequency summing to 1.
MLR has $2n$ parameters, so that for each variant $i$, we estimate its initial frequency $p_i$ as well as its fixed growth rate $f_i$.
Because growth rates are necessarily relative, we define an arbitrary "pivot" variant to compare to, where we fix $f_i=1$.

:::figure{#fig-time-vs-frequency-sarscov2 component=time-vs-frequency dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_frequency.png}
**Relative frequencies of SARS-CoV-2 clades through time.**
Points represent empirical frequencies of SARS-CoV-2 Nextstrain clades, while solid lines represent modeled frequencies from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

With variant frequency $x_i(t)$ and constant variant fitness $f_i$, we can describe the mean population fitness as a standard weighted sum $\bar{f}(t) = \sum_i x_i(t) \, f_i$.
The fitness flux [@mustonen2010fitness] of the population is then the rate of change of population fitness at a given time $\phi(t) = \Delta \bar{f}(t) / \Delta t$.



:::figure{#fig-time-vs-fitness-sarscov2 component=time-vs-fitness dataset=sarscov2_clades static=figures/sarscov2_clades_time_vs_fitness.png}
**Cumulative SARS-CoV-2 fitness flux.**
Empirical frequencies of SARS-CoV-2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

:::figure{#fig-time-vs-frequency-h3n2 component=time-vs-frequency dataset=h3n2_clades static=figures/h3n2_clades_time_vs_frequency.png}
**Relative frequencies of H3N2 clades through time.**
Points represent empirical frequencies of H3N2 Nextstrain clades, while solid lines represent modeled frequencies from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

:::figure{#fig-time-vs-fitness-h3n2 component=time-vs-fitness dataset=h3n2_clades static=figures/h3n2_clades_time_vs_fitness.png}
**Cumulative H3N2 fitness flux.**
Empirical frequencies of H3N2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

Fisher's fundamental theorem and the expected relationship of fitness variance and fitness flux.
Empirical investigation of this relationship in SARS-CoV-2 and H3N2.

:::figure{#fig-sarscov2-variance-flux src=figures/sarscov2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in SARS-CoV-2.**
Log fitness variance is compared to log fitness flux, where each dot represents a daily timepoint.
:::

:::figure{#fig-h3n2-variance-flux src=figures/h3n2_clades_fitness_variance_vs_flux.png}
**Fitness variance and fitness flux in H3N2.**
Log fitness variance is compared to log fitness flux, where each dot represents a daily timepoint.
:::

We shouldn't do a simple correlation of mutations against cumulative fitness flux due to phylogenetic non-independence.
Instead can rely on phylogenetic contrasts of parent and daughter lineages.
Pango lineages provide a convenient granular and hierarchical nomenclature well suited to this.

:::figure{#fig-delta-genome src=figures/sarscov2_lineage_delta_fitness_across_genome.png}
**Correlation of lineage-specific amino acid change to lineage-specific fitness change across regions of the SARS-CoV-2 genome.**
:::

:::figure{#fig-delta-time src=figures/sarscov2_lineage_delta_fitness_across_time.png}
**Correlation of lineage-specific amino acid change to lineage-specific fitness change over time.**
:::

:::figure{#fig-delta-evescape src=figures/sarscov2_lineage_delta_fitness_vs_evescape.png}
**Correlation of lineage-specific change in EvEscape score to lineage-specific fitness change.**
:::

## Conclusions

TBD

## Methods

TBD

## Acknowledgments

This research was supported in part by grant NSF PHY-2309135, the Gordon and Betty Moore Foundation grant no. 2919.02, and the Chan Zuckerberg Initiative DAF grant to the Kavli Institute for Theoretical Physics (KITP).
