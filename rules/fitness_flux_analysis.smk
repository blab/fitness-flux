"""
This part of the workflow runs the fitness-flux analysis: a Python port of
fitness-flux-analysis/fitness-flux.nb. For each analysis dataset it aggregates
the per-season MLR estimates into a scaffolded fitness scale, gathers empirical
frequencies and mean dates, characterizes the fitness wave (the "flux"), and
joins mutation counts to fitness. Outputs land in fitness-flux-analysis/results/
and are visualized by viz/fitness-flux.html.
"""

FITNESS_FLUX_ANALYSES = ["h3n2_clades", "sarscov2_clades", "sarscov2_lineages"]

FITNESS_FLUX_OUTPUTS = [
    "direct_fitness.tsv",
    "scaffolded_fitness.tsv",
    "frequencies.tsv",
    "mean_date.tsv",
    "flux_timeseries.tsv",
    "flux_summary.json",
    "mutation_fitness.tsv",
]


def _fitness_flux_season_inputs(wildcards):
    """The per-season MLR results that an analysis dataset aggregates."""
    members = [
        dataset
        for dataset in config["datasets"]
        if dataset.startswith(wildcards.analysis + "_")
    ]
    return expand("mlr-estimates/{member}/mlr_results.json", member=members)


rule fitness_flux_gather_fitness:
    input:
        mlr = _fitness_flux_season_inputs
    output:
        "fitness-flux-analysis/results/{analysis}_direct_fitness.tsv"
    log:
        "logs/fitness_flux/{analysis}_gather_fitness.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/gather_fitness.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            --output {output} 2>&1 | tee {log}
        """


rule fitness_flux_scaffold_fitness:
    input:
        mlr = _fitness_flux_season_inputs
    output:
        "fitness-flux-analysis/results/{analysis}_scaffolded_fitness.tsv"
    log:
        "logs/fitness_flux/{analysis}_scaffold_fitness.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/scaffold_fitness.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            --output {output} 2>&1 | tee {log}
        """


rule fitness_flux_gather_frequencies:
    input:
        mlr = _fitness_flux_season_inputs,
        scaffolded = "fitness-flux-analysis/results/{analysis}_scaffolded_fitness.tsv"
    output:
        frequencies = "fitness-flux-analysis/results/{analysis}_frequencies.tsv",
        mean_date = "fitness-flux-analysis/results/{analysis}_mean_date.tsv"
    log:
        "logs/fitness_flux/{analysis}_gather_frequencies.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/gather_frequencies.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            --scaffolded {input.scaffolded} \
            --frequencies-output {output.frequencies} \
            --mean-date-output {output.mean_date} 2>&1 | tee {log}
        """


rule fitness_flux_wave:
    input:
        scaffolded = "fitness-flux-analysis/results/{analysis}_scaffolded_fitness.tsv",
        frequencies = "fitness-flux-analysis/results/{analysis}_frequencies.tsv"
    output:
        timeseries = "fitness-flux-analysis/results/{analysis}_flux_timeseries.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_flux_summary.json"
    params:
        generation_time = config.get("analysis_generation_time", 3.2)
    log:
        "logs/fitness_flux/{analysis}_wave.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/fitness_wave.py \
            --dataset {wildcards.analysis} \
            --scaffolded {input.scaffolded} \
            --frequencies {input.frequencies} \
            --generation-time {params.generation_time} \
            --timeseries-output {output.timeseries} \
            --summary-output {output.summary} 2>&1 | tee {log}
        """


rule fitness_flux_mutation_fitness:
    input:
        mut_counts = "fitness-flux-analysis/source-data/{analysis}_mut_counts.tsv",
        scaffolded = "fitness-flux-analysis/results/{analysis}_scaffolded_fitness.tsv",
        mean_date = "fitness-flux-analysis/results/{analysis}_mean_date.tsv"
    output:
        "fitness-flux-analysis/results/{analysis}_mutation_fitness.tsv"
    log:
        "logs/fitness_flux/{analysis}_mutation_fitness.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/mutation_fitness.py \
            --mut-counts {input.mut_counts} \
            --scaffolded {input.scaffolded} \
            --mean-date {input.mean_date} \
            --output {output} 2>&1 | tee {log}
        """


rule all_fitness_flux:
    input:
        expand(
            "fitness-flux-analysis/results/{analysis}_{output}",
            analysis=FITNESS_FLUX_ANALYSES,
            output=FITNESS_FLUX_OUTPUTS,
        )
