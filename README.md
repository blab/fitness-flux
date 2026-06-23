# Comparing fitness dynamics across SARS-CoV-2, influenza H3 and influenza H1

## Workflow

Once metadata is provisioned locally, run the entire workflow with
```
nextstrain build . 
```

### Provision metadata

Input metadata is provisioned by the workflow itself (rule `provision_metadata`): it streams the per-virus Nextstrain metadata from S3, subsets to the columns the analysis uses, and writes `data/{virus}_subset_metadata.tsv.zst` (the `local_metadata` every dataset consumes). Provision both viruses with
```
nextstrain build . all_provision_metadata
```
or just run the full build (below) — provisioning runs automatically as an upstream dependency whenever the metadata files are missing. The source URLs and column sets live under `provision:` in `defaults/config.yaml`.

Notes:
- Requires `aws`, `zstd`, `xz`, and `tsv-select` on `PATH`.
- H3N2 reads a private bucket (`nextstrain-data-private`), so AWS credentials are required; SARS-CoV-2 uses the public `nextstrain-data` bucket.
- The files are provisioned once and not re-downloaded automatically; refresh with `nextstrain build . --forcerun provision_metadata all_provision_metadata` or by deleting `data/*_subset_metadata.tsv.zst`.

### Sequence counts

Data for the project consists of daily sequence counts of clades of SARS-CoV-2,
influenza H3 and influenza H1. Sequence counts are provisioned to the
`sequence-counts/` directory. From top-level directory run
```
nextstrain build . all_sequence_counts
```
to produce the sequence counts files
```
sequence-counts/sarscov2_clades_2020/prepared_seq_counts.tsv
```
Currently, clade counts are provisioned for just the USA.

### MLR estimates

Run MLR models using [evofr package](https://github.com/blab/evofr). Run the
model with
```
nextstrain build . all_mlr_estimates
```
to produce the MLR output JSON files
```
mlr-estimates/sarscov2_clades_2020/mlr_results.json
```

### Scaffolded fitnesses

Fitnesses within each timepoint are measured relative to an arbitrary pivot variant. The Mathematica notebook `fitness-flux-analysis/fitness-flux.nb` takes `mlr_results.json` across timepoints and combines into a single `scaffolded-fitness/sarscov2_clades_scaffolded_fitness.tsv`. This notebook needs to be run separately for virus `sarscov2` classification `clades`, virus `sarscov2` classification `lineages` and virus `h3n2` classification `clades`.
