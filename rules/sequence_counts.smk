"""
Summarize sequence counts from existing metadata for a particular dataset
"""

# rule subset_metadata:
#     input:
#         metadata = lambda w: config[w.dataset]["local_metadata"]
#     output:
#         subset_metadata = "sequence-counts/{dataset}/subset_metadata.tsv.zst"
#     params:
#         subset_columns = lambda w: ",".join(config[w.dataset]["subset_columns"]),
#     shell:
#         """
#         zstd -c -d {input.metadata} \
#             | tsv-select -H -f {params.subset_columns:q} \
#             | zstd -c > {output.subset_metadata}
#         """

rule summarize_clade_sequence_counts:
    input:
        subset_metadata = lambda w: config[w.dataset]["local_metadata"]
    output:
        clade_seq_counts = "sequence-counts/{dataset}/seq_counts.tsv"
    params:
        seq_count_options = lambda w: config[w.dataset]["seq_count_options"]
    shell:
        """
        ./scripts/summarize-clade-sequence-counts \
            --metadata {input.subset_metadata} \
            --output {output.clade_seq_counts} \
            {params.seq_count_options}
        """

def _get_prepare_data_option(wildcards, option_name):
    """
    Return the option for prepare data from the config based on the
    wildcards.dataset values.

    If the *option* exists as a key within config[wildcard.dataset]
    then return as "--{option-name} {option_value}". Or else return an empty string.
    """
    option_value = config.get(wildcards.dataset, {}) \
                         .get(option_name)

    if option_value is not None:
        # Change underscores of YAML keys to dashes for proper CLI option names
        option_name = option_name.replace('_', '-')
        return f'--{option_name} {option_value}'

    return ''

rule prepare_clade_data:
    """Preparing clade counts for analysis"""
    input:
        sequence_counts = "sequence-counts/{dataset}/seq_counts.tsv"
    output:
        sequence_counts = "sequence-counts/{dataset}/prepared_seq_counts.tsv"
    log:
        "logs/{dataset}/prepare_clade_data.txt"
    params:
        min_date = lambda wildcards: _get_prepare_data_option(wildcards, 'min_date'),
        max_date = lambda wildcards: _get_prepare_data_option(wildcards, 'max_date'),
        location_min_seq = lambda wildcards: _get_prepare_data_option(wildcards, 'location_min_seq'),
        location_min_seq_days = lambda wildcards: _get_prepare_data_option(wildcards, 'location_min_seq_days'),
        excluded_locations = lambda wildcards: _get_prepare_data_option(wildcards, 'excluded_locations'),
        clade_min_seq = lambda wildcards: _get_prepare_data_option(wildcards, 'clade_min_seq'),
        clade_min_seq_days = lambda wildcards: _get_prepare_data_option(wildcards, 'clade_min_seq_days'),
        force_include_clades = lambda wildcards: _get_prepare_data_option(wildcards, 'force_include_clades'),
        force_exclude_clades = lambda wildcards: _get_prepare_data_option(wildcards, 'force_exclude_clades')
    shell:
        """
        python ./scripts/prepare-data.py \
            --seq-counts {input.sequence_counts} \
            {params.min_date} \
            {params.max_date} \
            {params.location_min_seq} \
            {params.location_min_seq_days} \
            {params.excluded_locations} \
            {params.clade_min_seq} \
            {params.clade_min_seq_days} \
            {params.force_include_clades} \
            {params.force_exclude_clades} \
            --output-seq-counts {output.sequence_counts} 2>&1 | tee {log}
        """

rule collapse_sequence_counts:
    "Collapsing Pango lineages, based on sequence count threshold"
    input:
        sequence_counts = "sequence-counts/{dataset}/prepared_seq_counts.tsv",
    output:
        sequence_counts = "sequence-counts/{dataset}/collapsed_seq_counts.tsv"
    log:
        "logs/{dataset}/collapse_sequence_counts.txt"
    params:
        collapse_threshold = lambda wildcards: _get_prepare_data_option(wildcards, 'collapse_threshold'),
    shell:
        """
        python ./scripts/collapse-lineage-counts.py \
            --seq-counts {input.sequence_counts} \
            {params.collapse_threshold} \
            --output-seq-counts {output.sequence_counts} 2>&1 | tee {log}
        """

rule annotate_sequence_counts:
    input:
        sequence_counts = "sequence-counts/{dataset}/collapsed_seq_counts.tsv"
    output:
        sequence_counts = "sequence-counts/{dataset}/annotated_seq_counts.tsv"
    shell:
        """
        dataset_suffix=$(echo "{wildcards.dataset}" | awk -F'_' '{{print $NF}}')
        python scripts/annotate_sequence_counts.py \
            --input {input.sequence_counts} \
            --output {output.sequence_counts} \
            --dataset-suffix $dataset_suffix
        """

import os

# Rule to combine annotated sequence counts for each analysis
rule aggregate_sequence_counts:
    input:
        lambda wildcards: expand(
            "sequence-counts/{dataset}/annotated_seq_counts.tsv",
            dataset=[d for d in config["datasets"] if d.startswith(f"{wildcards.analysis}_")]
        )
    output:
        combined = "aggregated-counts/{analysis}/aggregated_sequence_counts.tsv"
    run:
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(output.combined), exist_ok=True)

        # Handle the case where no input files are present
        if len(input) == 0:
            with open(output.combined, 'w') as out_file:
                out_file.write("location\tvariant\tdate\tsequences\n")  # Example header
        else:
            shell("""
                (head -n 1 {input[0]} && tail -n +2 -q {input}) > {output.combined}
            """)
