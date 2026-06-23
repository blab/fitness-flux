# Comparing fitness dynamics across SARS-CoV-2, influenza H3 and influenza H1

## Provision metadata locally

```
mkdir data
cd data
```

For SARS-CoV-2
```
aws s3 cp s3://nextstrain-data/files/ncov/open/metadata.tsv.zst sarscov2_metadata.tsv.zst
zstd -c -d sarscov2_metadata.tsv.zst \
   | tsv-select -H -f strain,date,country,clade_nextstrain,Nextclade_pango,QC_overall_status \
   | zstd -c > data/sarscov2_subset_metadata.tsv.zst
```
This results in `data/sarscov2_subset_metadata.tsv.zst`

For H3N2
```
aws s3 cp s3://nextstrain-data-private/files/workflows/seasonal-flu/h3n2/metadata.tsv.xz h3n2_metadata.tsv.xz
xz -c -d h3n2_metadata.tsv.xz \
   | tsv-select -H -f strain,date,country,subclade_nextclade_ha,qc.overallStatus_ha \
   | zstd -c > data/h3n2_subset_metadata.tsv.zst
```
This results in `data/h3n2_subset_metadata.tsv.zst`

## Workflow

Once metadata is provisioned locally, run the entire workflow with
```
nextstrain build . all_mlr_estimates
```

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
