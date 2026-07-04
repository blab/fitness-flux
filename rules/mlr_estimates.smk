"""
This part of the workflow runs the model scripts.
"""

rule all_mlr_estimates:
    input:
        mlr_estimates = expand(
            "mlr-estimates/{virus}/mlr_results.json",
            virus=config["datasets"]
        )

def _get_models_option(wildcards, option_name):
    """
    Return the option for model from the config based on the
    wildcards.data_provenance, wildcards.variant_classification and the wildcards.geo_resolution values.

    If the *option* exists as a key within config['models'][wildcard.data_provenance][wildcard.variant_classification][wildcard.geo_resolution]
    then return as "--{option-name} {option_value}". Or else return an empty string.
    """
    option_value = config.get(wildcards.dataset, {}) \
                         .get(option_name)

    if option_value is not None:
        # Change underscores of YAML keys to dashes for proper CLI option names
        option_name = option_name.replace('_', '-')
        return f'--{option_name} {option_value}'

    return ''

def _generation_time_options(dataset):
    """Generation-time CLI flags from the virus-keyed config['generation_time'].

    SARS-CoV-2 is a dict carrying the per-variant pre/post-Omicron split (the split
    is applied per variant inside run-mlr-model.py); H3N2 is a single scalar. Also
    used by rules/fitness_flux_analysis.smk (lazy param lambdas, so include order
    does not matter).
    """
    virus = dataset.split("_")[0]
    gt = config["generation_time"][virus]
    if isinstance(gt, dict):
        classification = "lineages" if "lineages" in dataset else "clades"
        return (
            f"--generation-time {gt['post_omicron']}"
            f" --generation-time-pre-omicron {gt['pre_omicron']}"
            f" --variant-classification {classification}"
        )
    return f"--generation-time {gt}"

rule mlr_model:
    input:
        sequence_counts = "sequence-counts/{dataset}/collapsed_seq_counts.tsv"
    output:
        # Note this output is not used in the shell command because it is one of the many
        # files generated and output to the export path.
        # We are listing this specific file as the output file because it is the final
        # final output of the model script.
        results = "mlr-estimates/{dataset}/mlr_results.json"
    log:
        "logs/{dataset}/mlr_model.txt"
    params:
        model_config = config.get("mlr_config"),
        export_path = lambda w: f"mlr-estimates/{w.dataset}",
        pivot = lambda wildcards: _get_models_option(wildcards, 'pivot'),
        generation_time = lambda wildcards: _generation_time_options(wildcards.dataset),
        # Empirical weekly_raw_freq smoothing window (days), from config
        # (raw_freq_window). H3N2 sets 14; datasets without the key use the
        # run-mlr-model.py default of 7.
        raw_freq_window = lambda wildcards: _get_models_option(wildcards, 'raw_freq_window')
    resources:
        mem_mb=4000
    shell:
        """
        python -u ./scripts/run-mlr-model.py \
            --config {params.model_config} \
            --seq-path {input.sequence_counts} \
            --export-path {params.export_path} \
            {params.pivot} \
            {params.generation_time} \
            {params.raw_freq_window} \
            --data-name mlr 2>&1 | tee {log}
        """
