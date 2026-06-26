# Fitness flux in SARS-CoV-2 and influenza H3N2

Trevor Bedford<sup>1,2</sup>

1. Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA
2. Howard Hughes Medical Institute, Seattle, WA, USA

## Abstract

The tempo of viral adaptation is usually read indirectly from the composition of mutations, through measures such as dN/dS. Here we measure it directly from the dynamics of variant frequencies, where we use multinomial logistic regression to estimate a fitness for each co-circulating variant. We aggregate these estimates to derive the rate of change of mean population fitness, referred to as fitness flux. Tracing SARS-CoV-2 from 2020 to 2026 and comparing against seasonal influenza A/H3N2, we find that SARS-CoV-2 adapted rapidly with a 6.7-fold increase in fitness from 2020 to 2023, before slowing to a 2.2-fold increase from 2023 to 2026. Influenza H3N2 sustains a slower, steadier pace roughly threefold below recent SARS-CoV-2. In both, the rate of fitness gain closely tracks the variance in fitness, matching the 1:1 expectation of Fisher's fundamental theorem. Phylogenetic contrasts between parent and child lineages localize most fitness gain to spike, and within spike to the receptor-binding domain, where a simple count of spike S1 substitutions predicts lineage fitness about as well as deep-learning escape and protein-language-model scores. Measuring fitness directly thus offers a transparent, frequency-based alternative to mutational proxies for tracking and anticipating viral adaptation.

## Installation

Clone the repository and install Python dependencies:
```
pip install -r requirements.txt
```
`evofr` is needed to regenerate `mlr-estimates/` and the downstream analysis and visualization steps need  numpy/scipy/pandas.

Provisioning the input metadata additionally requires `aws`, `zstd`, `xz`, and `tsv-select` on `PATH`. H3N2 reads a private bucket (`nextstrain-data-private`), so AWS credentials are required; SARS-CoV-2 uses the public `nextstrain-data` bucket.

Workflow targets can be run either as `snakemake <target>` from this environment or as `nextstrain build . <target>` through the [Nextstrain CLI](https://docs.nextstrain.org/projects/cli/) runtime (`environment_nextstrain.sh`). The examples below use the latter.

## Workflow

Once metadata is provisioned locally, run the entire workflow with
```
nextstrain build .
```
The default target builds the full pipeline: sequence counts, variant relationships, MLR estimates, the fitness-flux analysis, and the lineage-deltas analysis. The individual stages can also be run on their own, as described below.

### Provision metadata

Input metadata is provisioned by the workflow itself (rule `provision_metadata`): it streams the per-virus Nextstrain metadata from S3, subsets to the columns the analysis uses, and writes `data/{virus}_subset_metadata.tsv.zst`. Provision both viruses with
```
nextstrain build . all_provision_metadata
```
or just run the full build — provisioning runs automatically as an upstream dependency whenever the metadata files are missing. The source URLs and column sets live under `provision:` in `defaults/config.yaml`. The files are provisioned once and not re-downloaded automatically; refresh with `nextstrain build . --forcerun provision_metadata all_provision_metadata` or by deleting `data/*_subset_metadata.tsv.zst`.

### Sequence counts

Daily clade and Pango-lineage sequence counts are provisioned to the `sequence-counts/` directory, with rare clades collapsed into "other" and rare lineages rolled up into their parents. Run
```
nextstrain build . all_sequence_counts all_variant_relationships
```
to produce, for each dataset, the collapsed counts and the lineage parent–child relationships
```
sequence-counts/sarscov2_clades_2020/collapsed_seq_counts.tsv
sequence-counts/sarscov2_lineages_2020/variant_relationships.tsv
```

### MLR estimates

Fit multinomial logistic regression (MLR) with the [evofr package](https://github.com/blab/evofr) on 1-year sliding windows for SARS-CoV-2 and 2-year windows for H3N2:
```
nextstrain build . all_mlr_estimates
```
to produce the MLR output JSON files
```
mlr-estimates/sarscov2_clades_2020/mlr_results.json
```

### Fitness-flux analysis

Within each window MLR measures fitness only relative to an arbitrary pivot. The `all_fitness_flux` target scaffolds the per-window estimates onto a single per-variant scale (rule `fitness_flux_scaffold_fitness`), then computes mean population fitness, fitness variance, and fitness flux through time:
```
nextstrain build . all_fitness_flux
```
This writes the analysis tables under `fitness-flux-analysis/results/` together with the per-figure data consumed by the interactive components in `viz/`
```
fitness-flux-analysis/results/sarscov2_clades_scaffolded_fitness.tsv
viz/time-vs-fitness/data/sarscov2_clades.json
```

### Lineage deltas

The `all_lineage_deltas` target counts amino-acid substitutions per Pango lineage (rule `lineage_mut_counts`), forms parent-to-child branches, and contrasts each branch's change in substitution count against its change in fitness — also comparing against the EvEscape and ESM-2 mutational-fitness predictors:
```
nextstrain build . all_lineage_deltas
```
producing the per-branch deltas and predictor correlations under `lineage-deltas-analysis/results/`
```
lineage-deltas-analysis/results/branch_deltas.tsv
lineage-deltas-analysis/results/predictor_correlations.tsv
```

### Cleaning

`nextstrain build . clean` removes everything the workflow generates (the `sequence-counts/`, `mlr-estimates/`, logs, analysis `results/`, and generated viz data) so the repo can be rerun from scratch, leaving the provisioned raw metadata in place. `nextstrain build . clean_analysis` removes only the fitness-flux and lineage-deltas intermediates, leaving the upstream sequence counts and MLR estimates intact.

## Organization

- `manuscript/` — paper source (`fitness-flux.md`, `fitness_flux.bib`, figures), built with `press`.
- `Snakefile`, `rules/` — workflow definition; each `rules/*.smk` file covers one pipeline stage.
- `scripts/` — sequence-count preparation and MLR scripts (`run-mlr-model.py`, `collapse-lineage-counts.py`, `prepare-pango-relationships.py`, …).
- `defaults/config.yaml` — datasets and per-dataset parameters (windows, collapse thresholds, model settings).
- `data/`, `sequence-counts/`, `mlr-estimates/` — provisioned metadata and generated workflow outputs (gitignored).
- `fitness-flux-analysis/` — scaffolding plus the fitness-flux and fitness-variance analysis, figures, and viz source data.
- `lineage-deltas-analysis/` — per-branch mutation-to-fitness contrasts and the EvEscape / ESM-2 predictor source data.
- `mutation-counts/` — per-lineage substitution counts by genome region.
- `inclusion-thresholds/` — collapse-threshold sweep and the rationale behind the chosen cutoffs.
- `viz/` — interactive figure components (HTML/JS) fed by the `viz/*/data/` tables the workflow generates.

## Citation

Manuscript in preparation; citation forthcoming.
