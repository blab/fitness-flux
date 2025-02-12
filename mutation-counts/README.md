# Mutation counts

Prepare mutation counts for SARS-CoV-2 lineages

## Prep SARS-CoV-2 lineages JSON

Download JSON
```
nextstrain remote download https://nextstrain.org/nextclade/sars-cov-2
```
This results in the local file `nextclade_sars-cov-2.json`.

Flatten this JSON:
```
python ../scripts/flatten_auspice_json.py --json nextclade_sars-cov-2.json --output nextclade_sars-cov-2_flat.json
```
to produce `nextclade_sars-cov-2_flat.json`

Run Mathematica notebook `mutation-counts-sarscov2-lineages.nb`. This produces `results/sarscov2_lineages_mut_counts.tsv`.

This file has been versioned to `lineage-deltas-analysis/source-data/sarscov2_lineages_mut_counts.tsv`. 