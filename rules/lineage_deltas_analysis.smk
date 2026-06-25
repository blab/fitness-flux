"""
This part of the workflow runs the lineage-deltas analysis: a Python port of
lineage-deltas-analysis/lineage-deltas.nb. It relates the change in fitness along
each SARS-CoV-2 parent->child lineage branch to the change in mutation counts and
in mutation-effect predictors (CoVFit, EvEscape, Bloom-lab DMS, ESM embeddings).
Outputs land in lineage-deltas-analysis/results/ and are visualized by
viz/lineage-deltas.html.
"""

LINEAGE_DELTAS_DATASET = "sarscov2_lineages"
LINEAGE_DELTAS_SOURCE = "lineage-deltas-analysis/source-data"

LINEAGE_DELTAS_SEASONS = [
    dataset
    for dataset in config["datasets"]
    if dataset.startswith(LINEAGE_DELTAS_DATASET + "_")
]


def _lineage_deltas_branch_inputs(wildcards):
    return {
        "mlr": expand(
            "mlr-estimates/{member}/mlr_results.json",
            member=LINEAGE_DELTAS_SEASONS,
        ),
        "relationships": expand(
            "sequence-counts/{member}/variant_relationships.tsv",
            member=LINEAGE_DELTAS_SEASONS,
        ),
    }


rule lineage_deltas_branches:
    input:
        unpack(_lineage_deltas_branch_inputs),
        mut_counts = "mutation-counts/results/sarscov2_lineages_mut_counts.tsv"
    output:
        "lineage-deltas-analysis/results/branch_deltas.tsv"
    log:
        "logs/lineage_deltas/branches.txt"
    shell:
        """
        python -u lineage-deltas-analysis/scripts/build_branch_deltas.py \
            --mlr-dir mlr-estimates \
            --seqcounts-dir sequence-counts \
            --mut-counts {input.mut_counts} \
            --output {output} 2>&1 | tee {log}
        """


rule lineage_deltas_predictors:
    input:
        unpack(_lineage_deltas_branch_inputs),
        covfit = f"{LINEAGE_DELTAS_SOURCE}/CoVFit_Predictions_fold_0_20231102v4.tsv",
        evescape = f"{LINEAGE_DELTAS_SOURCE}/COVID19_all.csv",
        dms = f"{LINEAGE_DELTAS_SOURCE}/clade_phenotypes.csv"
    output:
        "lineage-deltas-analysis/results/predictor_deltas.tsv"
    log:
        "logs/lineage_deltas/predictors.txt"
    shell:
        """
        python -u lineage-deltas-analysis/scripts/compute_predictor_deltas.py \
            --mlr-dir mlr-estimates \
            --seqcounts-dir sequence-counts \
            --covfit {input.covfit} \
            --evescape {input.evescape} \
            --dms {input.dms} \
            --output {output} 2>&1 | tee {log}
        """


rule lineage_deltas_esm:
    input:
        unpack(_lineage_deltas_branch_inputs),
        embeddings = expand(
            f"{LINEAGE_DELTAS_SOURCE}/embeddings_{{model}}.tsv.xz",
            model=["650M_pretrained", "650M_fine_tuned", "3B_pretrained", "3B_fine_tuned"],
        )
    output:
        "lineage-deltas-analysis/results/esm_deltas.tsv"
    log:
        "logs/lineage_deltas/esm.txt"
    shell:
        """
        python -u lineage-deltas-analysis/scripts/compute_esm_deltas.py \
            --mlr-dir mlr-estimates \
            --seqcounts-dir sequence-counts \
            --source-data {LINEAGE_DELTAS_SOURCE} \
            --output {output} 2>&1 | tee {log}
        """


rule lineage_deltas_linear_models:
    input:
        branch_deltas = "lineage-deltas-analysis/results/branch_deltas.tsv"
    output:
        coefficients = "lineage-deltas-analysis/results/linear_model_coefficients.tsv",
        slope = "lineage-deltas-analysis/results/slope_through_time.tsv"
    log:
        "logs/lineage_deltas/linear_models.txt"
    shell:
        """
        python -u lineage-deltas-analysis/scripts/linear_models.py \
            --branch-deltas {input.branch_deltas} \
            --coefficients-output {output.coefficients} \
            --slope-output {output.slope} 2>&1 | tee {log}
        """


rule lineage_deltas_correlations:
    input:
        branch_deltas = "lineage-deltas-analysis/results/branch_deltas.tsv",
        predictor_deltas = "lineage-deltas-analysis/results/predictor_deltas.tsv",
        esm_deltas = "lineage-deltas-analysis/results/esm_deltas.tsv"
    output:
        "lineage-deltas-analysis/results/predictor_correlations.tsv"
    log:
        "logs/lineage_deltas/correlations.txt"
    shell:
        """
        python -u lineage-deltas-analysis/scripts/predictor_correlations.py \
            --branch-deltas {input.branch_deltas} \
            --predictor-deltas {input.predictor_deltas} \
            --esm-deltas {input.esm_deltas} \
            --output {output} 2>&1 | tee {log}
        """


rule viz_lineage_deltas_data:
    input:
        branch_deltas = "lineage-deltas-analysis/results/branch_deltas.tsv",
        predictor_deltas = "lineage-deltas-analysis/results/predictor_deltas.tsv",
        esm_deltas = "lineage-deltas-analysis/results/esm_deltas.tsv"
    output:
        deltas = "viz/lineage-deltas/data/sarscov2_lineages.json",
        trends = "viz/lineage-delta-trends/data/sarscov2_lineages.json",
        histograms = "viz/lineage-delta-histograms/data/sarscov2_lineages.json"
    log:
        "logs/lineage_deltas/viz_lineage_deltas_data.txt"
    shell:
        """
        python -u lineage-deltas-analysis/scripts/viz_lineage_deltas_data.py \
            --branch-deltas {input.branch_deltas} \
            --predictor-deltas {input.predictor_deltas} \
            --esm-deltas {input.esm_deltas} \
            --output {output.deltas} 2>&1 | tee {log}
        cp {output.deltas} {output.trends}
        cp {output.deltas} {output.histograms}
        """


rule all_lineage_deltas:
    input:
        "lineage-deltas-analysis/results/branch_deltas.tsv",
        "lineage-deltas-analysis/results/predictor_deltas.tsv",
        "lineage-deltas-analysis/results/esm_deltas.tsv",
        "lineage-deltas-analysis/results/linear_model_coefficients.tsv",
        "lineage-deltas-analysis/results/slope_through_time.tsv",
        "lineage-deltas-analysis/results/predictor_correlations.tsv",
        "viz/lineage-deltas/data/sarscov2_lineages.json",
        "viz/lineage-delta-trends/data/sarscov2_lineages.json"
