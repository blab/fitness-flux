"""
Forecast-accuracy analysis (Abousamra et al. Fig 2B adaptation): mean absolute
error of MLR vs naive clade-frequency forecasts as a function of forecast lead
time, for all four virus lineages. Downstream of the already-computed per-window
MLR fits (mlr-estimates/); it does not re-run fitting. Per-lineage results land
in fitness-flux-analysis/results/ and are consolidated into a single four-panel
figure (viz/forecast-accuracy/, data/all.json).

Included after rules/mlr_estimates.smk and rules/fitness_flux_analysis.smk so the
_generation_time_options helper and the {analysis} wildcard constraint are in
scope.
"""

FORECAST_ANALYSES = ["sarscov2_clades", "h3n2_clades", "h1n1pdm_clades", "vic_clades"]

# Panel order + display labels for the combined figure.
FORECAST_PANELS = [
    ("sarscov2_clades", "SARS-CoV-2"),
    ("h3n2_clades", "H3N2"),
    ("h1n1pdm_clades", "H1N1pdm"),
    ("vic_clades", "B/Victoria"),
]

# All clade datasets now use 1-year windows sliding 3 months. To keep a 6-month
# forecast horizon, pair each window with the one two quarters (~180 days) ahead as
# the truth source (the forecast_accuracy.py gap tolerance rejects the adjacent
# quarter at 91 days). Uniform across lineages.
FORECAST_SLIDE_DAYS = 180


def _forecast_slide(analysis):
    return FORECAST_SLIDE_DAYS


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
        # Generation time (helper in mlr_estimates.smk): SARS-CoV-2 expands to the
        # per-variant pre/post-Omicron flags, flu to a single tau. Matches the tau
        # used at fit time so delta = ln(ga)/tau is recovered exactly.
        generation_time = lambda wildcards: _generation_time_options(wildcards.analysis),
        slide = lambda wildcards: _forecast_slide(wildcards.analysis)
    log:
        "logs/fitness_flux/{analysis}_forecast_accuracy.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/forecast_accuracy.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            --slide-days {params.slide} \
            {params.generation_time} \
            --detail-output {output.detail} \
            --summary-output {output.summary} 2>&1 | tee {log}
        """


def _forecast_panel_args(wildcards):
    """`--panel KEY LABEL SUMMARY DETAIL` groups, one per lineage, in panel order."""
    parts = []
    for key, label in FORECAST_PANELS:
        parts.append(
            f'--panel {key} "{label}" '
            f'fitness-flux-analysis/results/{key}_forecast_accuracy_summary.json '
            f'fitness-flux-analysis/results/{key}_forecast_accuracy.tsv'
        )
    return " ".join(parts)


rule viz_forecast_accuracy_data:
    """Consolidate every lineage's forecast results into the combined four-panel
    data/all.json (one panel per lineage) plus the component meta.json."""
    input:
        summaries = expand(
            "fitness-flux-analysis/results/{analysis}_forecast_accuracy_summary.json",
            analysis=FORECAST_ANALYSES,
        ),
        details = expand(
            "fitness-flux-analysis/results/{analysis}_forecast_accuracy.tsv",
            analysis=FORECAST_ANALYSES,
        )
    output:
        data = "viz/forecast-accuracy/data/all.json",
        meta = "viz/forecast-accuracy/meta.json"
    params:
        panels = _forecast_panel_args
    log:
        "logs/fitness_flux/viz_forecast_accuracy.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_forecast_accuracy_data.py \
            {params.panels} \
            --output {output.data} \
            --meta-output {output.meta} 2>&1 | tee {log}
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
        "viz/forecast-accuracy/data/all.json",
        "viz/forecast-accuracy/meta.json"
