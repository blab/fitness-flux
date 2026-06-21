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

Description of MLR model and cumulative fitness flux.
Rate of evolution of SARS-CoV-2 in terms of fitness flux and in terms of S1 mutations.
Rate of evolution of H3N2 in terms of fitness flux and in terms of S1 mutations.

:::figure{#fig-sarscov2-frequencies src=figures/sarscov2_clades_frequencies.png}
**Relative frequencies of SARS-CoV-2 clades through time.**
Points represent empirical frequencies of SARS-CoV-2 Nextstrain clades, while solid lines represent modeled frequencies from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

:::figure{#fig-sarscov2-fitness-mutations src=figures/sarscov2_clades_fitnesses_mutations.png}
**Cumulative SARS-CoV-2 fitness flux and spike S1 mutations.**
(A) Empirical frequencies of SARS-CoV-2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
(B) As before, empirical frequencies of SARS-CoV-2 clades are represented by vertical thickness, though placement on the y-axis represents cumulative median spike S1 mutations in viruses belonging to each clade.
(C) Cumulative spike S1 mutations plotted against cumulative fitness flux across SARS-CoV-2 clades.
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

:::figure{#fig-h3n2-frequencies src=figures/h3n2_clades_frequencies.png}
**Relative frequencies of H3N2 clades through time.**
Points represent empirical frequencies of H3N2 Nextstrain clades, while solid lines represent modeled frequencies from Multinomial Logistic Regression (MLR).
All data is taken from the USA.
The MLR analysis assumes that the fitness of each clade is constant through time.
:::

:::figure{#fig-h3n2-fitness-mutations src=figures/h3n2_clades_fitnesses_mutations.png}
**Cumulative H3N2 fitness flux and spike S1 mutations.**
(A) Empirical frequencies of H3N2 Nextstrain clades are represented by vertical thickness and placement on the y-axis represents cumulative fitness flux estimated from Multinomial Logistic Regression (MLR).
(B) As before, empirical frequencies of H3N2 clades are represented by vertical thickness, though placement on the y-axis represents cumulative median spike S1 mutations in viruses belonging to each clade.
(C) Cumulative spike S1 mutations plotted against cumulative fitness flux across H3N2 clades.
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
