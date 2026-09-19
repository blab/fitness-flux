"""
Forecast-accuracy analysis (Abousamra et al. Fig 2B adaptation): mean absolute
error of MLR vs naive clade-frequency forecasts as a function of forecast lead
time. Downstream of the already-computed per-window MLR fits (mlr-estimates/);
it does not re-run fitting. Outputs land in fitness-flux-analysis/results/ and
are visualized by viz/forecast-accuracy/.

Included after rules/mlr_estimates.smk and rules/fitness_flux_analysis.smk so the
_generation_time_options helper and the {analysis} wildcard constraint are in
scope.
"""

# SARS-CoV-2 clades only for now. The analysis is dataset-agnostic (flu also
# slides on fixed windows), so extend this list to enable flu later.
FORECAST_ANALYSES = ["sarscov2_clades"]


def _forecast_season_inputs(wildcards):
    """The per-window MLR results an analysis dataset aggregates (fan-in)."""
    members = [
        dataset
        for dataset in config["datasets"]
        if dataset.startswith(wildcards.analysis + "_")
    ]
    return expand("mlr-estimates/{member}/mlr_results.json", member=members)


rule forecast_accuracy:
    input:
        mlr = _forecast_season_inputs
    output:
        detail = "fitness-flux-analysis/results/{analysis}_forecast_accuracy.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_forecast_accuracy_summary.json"
    params:
        # Per-variant pre/post-Omicron generation time (helper in mlr_estimates.smk),
        # matching the tau used at fit time so delta = ln(ga)/tau is recovered exactly.
        generation_time = lambda wildcards: _generation_time_options(wildcards.analysis)
    log:
        "logs/fitness_flux/{analysis}_forecast_accuracy.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/forecast_accuracy.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            {params.generation_time} \
            --detail-output {output.detail} \
            --summary-output {output.summary} 2>&1 | tee {log}
        """


rule viz_forecast_accuracy_data:
    """Build the forecast-accuracy component's data.json (per dataset): the
    aggregate MAE(lead) curve plus faint per-window-pair curves, errors in %."""
    input:
        detail = "fitness-flux-analysis/results/{analysis}_forecast_accuracy.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_forecast_accuracy_summary.json"
    output:
        "viz/forecast-accuracy/data/{analysis}.json"
    log:
        "logs/fitness_flux/{analysis}_viz_forecast_accuracy.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_forecast_accuracy_data.py \
            --detail {input.detail} \
            --summary {input.summary} \
            --output {output} 2>&1 | tee {log}
        """


rule viz_forecast_accuracy_meta:
    """Emit the forecast-accuracy component's meta.json (only the datasets it
    carries; currently just sarscov2_clades)."""
    output:
        "viz/forecast-accuracy/meta.json"
    params:
        include = ",".join(FORECAST_ANALYSES)
    log:
        "logs/fitness_flux/viz_forecast_accuracy_meta.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_meta.py \
            --include {params.include} \
            --output {output} 2>&1 | tee {log}
        """


rule all_forecast_analysis:
    input:
        expand(
            "fitness-flux-analysis/results/{analysis}_forecast_accuracy.tsv",
            analysis=FORECAST_ANALYSES,
        ),
        expand(
            "fitness-flux-analysis/results/{analysis}_forecast_accuracy_summary.json",
            analysis=FORECAST_ANALYSES,
        ),
        expand(
            "viz/forecast-accuracy/data/{analysis}.json",
            analysis=FORECAST_ANALYSES,
        ),
        "viz/forecast-accuracy/meta.json"
