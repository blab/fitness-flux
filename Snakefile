configfile: "defaults/config.yaml"

rule all_sequence_counts:
    input:
        sequence_counts = expand(
            "sequence-counts/{virus}/collapsed_seq_counts.tsv",
            virus=config["datasets"]
        )

rule all_variant_relationships:
    input:
        variant_relationships = expand(
            "sequence-counts/{virus}/variant_relationships.tsv",
            virus=config["datasets"]
        )

rule all_mlr_estimates:
    input:
        mlr_estimates = expand(
            "mlr-estimates/{virus}/mlr_results.json",
            virus=config["datasets"]
        )

rule all:
    input:
        sequence_counts = expand(
            "sequence-counts/{virus}/collapsed_seq_counts.tsv",
            virus=config["datasets"]
        ),
        variant_relationships = expand(
            "sequence-counts/{virus}/variant_relationships.tsv",
            virus=config["datasets"]
        ),
        mlr_estimates = expand(
            "mlr-estimates/{virus}/mlr_results.json",
            virus=config["datasets"]
        )

include: "rules/sequence_counts.smk"
include: "rules/mlr_estimates.smk"
