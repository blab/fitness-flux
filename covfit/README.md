# Analysis of CoVFit fitness estimates

CoVFit is described in [Ito et al. A protein language model for exploring viral fitness landscapes. Nature Comms](https://www.nature.com/articles/s41467-025-59422-w)
with associated GitHub repo [github.com/TheSatoLab/CoVFit](https://github.com/TheSatoLab/CoVFit).

## Download CoVFit CLI

Start by downloading version `20231102v4` from https://zenodo.org/records/14614868. I believe this is a model trained on data up to 2023-11-02.

After decompressing this gives the directory `CoVFit_CLI` move this directory within this repo directly under `covfit/`

## Provision lineage spike sequences

Download JSON
```
nextstrain remote download https://nextstrain.org/nextclade/sars-cov-2
```

Extract sequences

```
python ../scripts/alignment.py --tree nextclade_sars-cov-2.json --root nextclade_sars-cov-2_root-sequence.json --output results/lineages.fasta --gene S --tips-only
```

# Run CoVFit estimates on these sequences

```
python CoVFit_CLI/run_covfit.py --input results/lineages.fasta
```

This produces the file `CoVFit_Predictions_fold_0.tsv` take this file and rename to `CoVFit_Predictions_fold_0_20231102v4.tsv` to include model name and version to `lineage-deltas-analysis/source-data/`.
