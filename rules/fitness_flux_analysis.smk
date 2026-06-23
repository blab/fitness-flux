"""
This part of the workflow runs the fitness-flux analysis: a Python port of
fitness-flux-analysis/fitness-flux.nb. For each analysis dataset it aggregates
the per-season MLR estimates into a scaffolded fitness scale, gathers empirical
frequencies and mean dates, characterizes the fitness wave (the "flux"), and
joins mutation counts to fitness. Outputs land in fitness-flux-analysis/results/
and are visualized by viz/fitness-flux.html.
"""

FITNESS_FLUX_ANALYSES = ["h3n2_clades", "sarscov2_clades", "sarscov2_lineages"]

# Pin {analysis} to the dataset names so the greedy wildcard can't, e.g., claim
# "{analysis}_frequencies.tsv" for "..._seasonal_frequencies.tsv".
wildcard_constraints:
    analysis = "|".join(FITNESS_FLUX_ANALYSES)

FITNESS_FLUX_OUTPUTS = [
    "direct_fitness.tsv",
    "scaffolded_fitness.tsv",
    "frequencies.tsv",
    "mean_date.tsv",
    "flux_timeseries.tsv",
    "flux_summary.json",
    "mutation_fitness.tsv",
    "colors.tsv",
    "seasonal_frequencies.tsv",
]

# The Nextstrain tree JSON supplies the sarscov2_clades clade colors; the other
# datasets derive colors from a Rainbow gradient and need no external input.
NCOV_TREE_JSON = "fitness-flux-analysis/ncov_global_all-time_6k_2026-03-02.json"


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


def _fitness_flux_colors_inputs(wildcards):
    inputs = {
        "scaffolded": f"fitness-flux-analysis/results/{wildcards.analysis}_scaffolded_fitness.tsv",
        "mean_date": f"fitness-flux-analysis/results/{wildcards.analysis}_mean_date.tsv",
        "frequencies": f"fitness-flux-analysis/results/{wildcards.analysis}_frequencies.tsv",
    }
    if wildcards.analysis == "sarscov2_clades":
        inputs["ncov"] = NCOV_TREE_JSON
    return inputs


rule fitness_flux_colors:
    input:
        unpack(_fitness_flux_colors_inputs)
    output:
        "fitness-flux-analysis/results/{analysis}_colors.tsv"
    params:
        ncov_flag = lambda w: f"--ncov-json {NCOV_TREE_JSON}" if w.analysis == "sarscov2_clades" else ""
    log:
        "logs/fitness_flux/{analysis}_colors.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/colors.py \
            --dataset {wildcards.analysis} \
            --scaffolded {input.scaffolded} \
            --mean-date {input.mean_date} \
            --frequencies {input.frequencies} \
            {params.ncov_flag} \
            --output {output} 2>&1 | tee {log}
        """


rule fitness_flux_seasonal_frequencies:
    input:
        mlr = _fitness_flux_season_inputs
    output:
        "fitness-flux-analysis/results/{analysis}_seasonal_frequencies.tsv"
    log:
        "logs/fitness_flux/{analysis}_seasonal_frequencies.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/seasonal_frequencies.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            --output {output} 2>&1 | tee {log}
        """


rule viz_fitness_flux_data:
    """Build the data.json (per dataset) for the time-vs-fitness AND
    frequency-vs-fitness components: the frequency/fitness join plus the embedded
    variant color table. The two components consume identical data (one plots it
    against time, the other against frequency), so build once and copy."""
    input:
        frequencies = "fitness-flux-analysis/results/{analysis}_frequencies.tsv",
        scaffolded = "fitness-flux-analysis/results/{analysis}_scaffolded_fitness.tsv",
        colors = "fitness-flux-analysis/results/{analysis}_colors.tsv"
    output:
        fitness = "viz/time-vs-fitness/data/{analysis}.json",
        freq_vs_fitness = "viz/frequency-vs-fitness/data/{analysis}.json"
    log:
        "logs/fitness_flux/{analysis}_viz_fitness_flux.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_fitness_flux_data.py \
            --frequencies {input.frequencies} \
            --scaffolded {input.scaffolded} \
            --colors {input.colors} \
            --output {output.fitness} 2>&1 | tee {log}
        cp {output.fitness} {output.freq_vs_fitness}
        """


rule viz_frequency_panels_data:
    """Build the time-vs-frequency component's data.json (per dataset): the
    per-season empirical/modeled frequencies plus the embedded color table."""
    input:
        seasonal = "fitness-flux-analysis/results/{analysis}_seasonal_frequencies.tsv",
        colors = "fitness-flux-analysis/results/{analysis}_colors.tsv"
    output:
        "viz/time-vs-frequency/data/{analysis}.json"
    log:
        "logs/fitness_flux/{analysis}_viz_frequency_panels.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_frequency_panels_data.py \
            --seasonal {input.seasonal} \
            --colors {input.colors} \
            --output {output} 2>&1 | tee {log}
        """


rule viz_variance_flux_data:
    """Build the variance-vs-flux component's data.json (per dataset): the daily
    fitness-wave timeseries (variance + velocity) plus the variance-vs-velocity
    regression, for the four-panel fitness-wave figure."""
    input:
        timeseries = "fitness-flux-analysis/results/{analysis}_flux_timeseries.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_flux_summary.json"
    output:
        "viz/variance-vs-flux/data/{analysis}.json"
    log:
        "logs/fitness_flux/{analysis}_viz_variance_flux.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_variance_flux_data.py \
            --timeseries {input.timeseries} \
            --summary {input.summary} \
            --output {output} 2>&1 | tee {log}
        """


rule viz_meta:
    """Emit each component's meta.json manifest (dataset ids + labels + default)
    for the dashboard selector and dev harness. Content is shared, so write both."""
    output:
        fitness_flux = "viz/time-vs-fitness/meta.json",
        frequency_panels = "viz/time-vs-frequency/meta.json",
        frequency_vs_fitness = "viz/frequency-vs-fitness/meta.json",
        variance_vs_flux = "viz/variance-vs-flux/meta.json"
    log:
        "logs/fitness_flux/viz_meta.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_meta.py --output {output.fitness_flux} 2>&1 | tee {log}
        python -u fitness-flux-analysis/scripts/viz_meta.py --output {output.frequency_panels} 2>&1 | tee -a {log}
        python -u fitness-flux-analysis/scripts/viz_meta.py --output {output.frequency_vs_fitness} 2>&1 | tee -a {log}
        python -u fitness-flux-analysis/scripts/viz_meta.py --output {output.variance_vs_flux} 2>&1 | tee -a {log}
        """


rule all_fitness_flux:
    input:
        expand(
            "fitness-flux-analysis/results/{analysis}_{output}",
            analysis=FITNESS_FLUX_ANALYSES,
            output=FITNESS_FLUX_OUTPUTS,
        ),
        expand(
            "viz/time-vs-fitness/data/{analysis}.json",
            analysis=FITNESS_FLUX_ANALYSES,
        ),
        expand(
            "viz/time-vs-frequency/data/{analysis}.json",
            analysis=FITNESS_FLUX_ANALYSES,
        ),
        expand(
            "viz/frequency-vs-fitness/data/{analysis}.json",
            analysis=FITNESS_FLUX_ANALYSES,
        ),
        expand(
            "viz/variance-vs-flux/data/{analysis}.json",
            analysis=FITNESS_FLUX_ANALYSES,
        ),
        "viz/time-vs-fitness/meta.json",
        "viz/time-vs-frequency/meta.json",
        "viz/frequency-vs-fitness/meta.json",
        "viz/variance-vs-flux/meta.json"
