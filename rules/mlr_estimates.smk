"""
This part of the workflow runs the model scripts.
"""

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
        generation_time = lambda wildcards: _get_models_option(wildcards, 'generation_time')
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
            --data-name mlr 2>&1 | tee {log}
        """
