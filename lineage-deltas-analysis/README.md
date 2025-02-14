# Analysis of changes in mutations and fitness along Pango branches

## Prep

Running the analysis notebook `lineage-deltas.nb` requires timepoint-specific MLR estimates populated as:
```
mlr-estimates/sarscov2_lineages_2020/mlr_results.json
mlr-estimates/sarscov2_lineages_2020-21/mlr_results.json
...
```
as well as timepoint-specific Pango lineage relationships populated as:
```
sequence-counts/sarscov2_lineages_2020/variant_relationships.tsv
sequence-counts/sarscov2_lineages_2020-21/variant_relationships.tsv
...
```

Additionally, source data has been versioned to `lineage-deltas-analysis/source-data/` containing:
- `sarscov2_lineages_mut_counts.tsv` versioned after running notebook in `mutation-counts/`
- `COVID19_all.csv` of EvEscape scores downloaded via "All Strain Data" at https://evescape.org/data
- `embeddings_fine_tuned.tsv` and `log_likelihoods_fine_tuned.tsv` from running ESM-2 on SARS-CoV-2 spike sequences (experimental)
