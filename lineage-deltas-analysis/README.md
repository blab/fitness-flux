# Analysis of changes in mutations and fitness along Pango branches

This analysis relates the change in fitness along each SARS-CoV-2
parent&rarr;child lineage branch to the mutations acquired and to mutation-effect
predictors (CoVFit, EvEscape, Bloom-lab deep mutational scanning, and ESM
protein-language-model embeddings).

The calculations were originally done in the Mathematica notebook
`lineage-deltas.nb` (retained for reference). They are now reimplemented as Python
scripts in `scripts/`, wired into Snakemake via
`rules/lineage_deltas_analysis.smk`, writing TSV to `results/`. The narrative is
presented interactively by `viz/lineage-deltas.html`.

## Inputs

Timepoint-specific MLR estimates and Pango lineage relationships from the main
pipeline:
```
mlr-estimates/sarscov2_lineages_2020/mlr_results.json
mlr-estimates/sarscov2_lineages_2020-21/mlr_results.json
...
sequence-counts/sarscov2_lineages_2020/variant_relationships.tsv
sequence-counts/sarscov2_lineages_2020-21/variant_relationships.tsv
...
```
A branch's fitness change uses that season's `log(growth advantage)` for the
child minus the parent.

Source data versioned in `source-data/`:
- `sarscov2_lineages_mut_counts.tsv` — per-region mutation counts (from the notebook in `mutation-counts/`)
- `COVID19_all.csv` — EvEscape scores, "All Strain Data" at https://evescape.org/data
- `CoVFit_Predictions_fold_0_20231102v4.tsv` — CoVFit predicted fitness
- `clade_phenotypes.csv` — Bloom-lab DMS phenotypes, https://github.com/jbloomlab/SARS2-spike-predictor-phenos
- `embeddings_{650M,3B}_{pretrained,fine_tuned}.tsv.xz` — ESM-2 embeddings, https://github.com/blab/embedded-pathways

## Pipeline

From the repository root:
```
pip install -r requirements.txt
snakemake all_lineage_deltas
```

The rules run these scripts:

| Script | Output(s) in `results/` | Purpose |
| --- | --- | --- |
| `build_branch_deltas.py` | `branch_deltas.tsv` | per-branch change in per-region mutation counts and log fitness |
| `compute_predictor_deltas.py` | `predictor_deltas.tsv` | per-branch change in CoVFit / EvEscape / DMS predictors (long format) |
| `compute_esm_deltas.py` | `esm_deltas.tsv` | Euclidean distance between parent and child ESM embeddings |
| `linear_models.py` | `linear_model_coefficients.tsv`, `linear_model_predictions.tsv`, `slope_through_time.tsv` | multiple regression of fitness change on non-overlapping regions (RBD, S1-excl-RBD, ORF1ab, accessory) with estimate/SE/t/p; per-branch fitted values; per-season slopes |
| `predictor_correlations.py` | `predictor_correlations.tsv` | Pearson / Spearman / R² of each predictor delta vs. fitness change |

## Visualization

`viz/lineage-deltas.html` (a self-contained Observable Plot / d3 page served as
part of the static GitHub Pages site) fetches these `results/` tables at runtime
and renders the branch mutation–fitness scatter, slope through time, the linear
model coefficients, and the predictor correlations. The PNGs under `figures/`
remain from the original notebook.
