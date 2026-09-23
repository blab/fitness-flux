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


# ---------------------------------------------------------------------------
# Similarity-aware scoring (H3N2)
#
# The label-based score above treats clade labels as exchangeable, so a forecast
# that lands on a clade's parent is charged exactly as much as one that lands on
# a distant clade. These rules score the SAME forecasts in HA1 mutation-profile
# space, where near-misses cost less, using one global clade similarity map built
# from the reconstructed HA1 sequence at each clade's MRCA node.
#
# H3N2 only: the map needs a per-dataset auspice tree with a gene annotation and
# clade labels matching the analysis vocabulary, and the vaccine-strain-selection
# framing this serves is influenza-specific.
# ---------------------------------------------------------------------------

SIMILARITY_ANALYSES = [a for a in FORECAST_ANALYSES if a in config.get("clade_distance", {}).get("trees", {})]

# Both the existing 6-month pairing and a 1-year pairing. Clade labels drift more
# the further out you forecast, so the longer horizon is where the label-based
# score loses the most information — and it is the vaccine-selection timeline.
SIMILARITY_HORIZONS = [180, 365]


def _clade_tree(analysis, field):
    return config["clade_distance"]["trees"][analysis][field]


rule download_clade_tree:
    """Fetch the seasonal-flu auspice tree that defines the clade similarity map.

    The published dataset embeds its own root sequence, so no sidecar
    root-sequence file is needed. Not re-downloaded automatically; refresh with
    `--forcerun download_clade_tree`.
    """
    output:
        tree = "data/{analysis}_clade_tree.json"
    params:
        url = lambda wildcards: _clade_tree(wildcards.analysis, "url")
    log:
        "logs/fitness_flux/{analysis}_download_clade_tree.txt"
    shell:
        """
        curl -sSL --compressed --fail --max-time 600 {params.url:q} -o {output.tree} 2>&1 | tee {log}
        """


rule clade_distances:
    """Reconstruct each clade's MRCA gene sequence and build the global map.

    Emits the per-clade sequences, the clade x (position, amino acid) indicator
    table, the pairwise amino-acid distance matrix, and a summary carrying the
    diagnostics the figure and Methods quote.
    """
    input:
        tree = "data/{analysis}_clade_tree.json",
        mlr = _forecast_season_inputs
    output:
        sequences = "fitness-flux-analysis/results/{analysis}_ha1_mrca.tsv",
        distances = "fitness-flux-analysis/results/{analysis}_distance_matrix.tsv",
        profile = "fitness-flux-analysis/results/{analysis}_profile_table.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_distance_summary.json"
    params:
        gene = lambda wildcards: _clade_tree(wildcards.analysis, "gene"),
        clade_key = lambda wildcards: _clade_tree(wildcards.analysis, "clade_key")
    log:
        "logs/fitness_flux/{analysis}_clade_distances.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/clade_distances.py \
            --tree {input.tree} \
            --gene {params.gene} \
            --clade-key {params.clade_key} \
            --mlr-dir mlr-estimates \
            --dataset {wildcards.analysis} \
            --sequences-output {output.sequences} \
            --distances-output {output.distances} \
            --profile-output {output.profile} \
            --summary-output {output.summary} 2>&1 | tee {log}
        """


rule forecast_similarity_selftest:
    """Assert the scoring transform is correct before any number is produced.

    Checks that a large eps reproduces the label-based score, that two pure
    clades sit exactly their amino-acid distance apart, and that shifting
    forecast mass to a more distant clade always costs more. A sentinel file
    makes the check a hard dependency of the scoring rule rather than a habit.
    """
    input:
        profile = "fitness-flux-analysis/results/{analysis}_profile_table.tsv",
        distances = "fitness-flux-analysis/results/{analysis}_distance_matrix.tsv"
    output:
        sentinel = touch("fitness-flux-analysis/results/{analysis}_similarity_selftest.ok")
    log:
        "logs/fitness_flux/{analysis}_similarity_selftest.txt"
    shell:
        """
        for norm in l1 l2; do
            python -u fitness-flux-analysis/scripts/forecast_similarity.py \
                --dataset {wildcards.analysis} \
                --profile-table {input.profile} \
                --distance-matrix {input.distances} \
                --norm $norm --self-test \
                --detail-output /dev/null --summary-output /dev/null
        done 2>&1 | tee {log}
        """


rule forecast_similarity:
    input:
        profile = "fitness-flux-analysis/results/{analysis}_profile_table.tsv",
        distances = "fitness-flux-analysis/results/{analysis}_distance_matrix.tsv",
        selftest = "fitness-flux-analysis/results/{analysis}_similarity_selftest.ok",
        mlr = _forecast_season_inputs
    output:
        detail = "fitness-flux-analysis/results/{analysis}_similarity_{horizon}d.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_similarity_{horizon}d_summary.json"
    wildcard_constraints:
        horizon = "|".join(str(h) for h in SIMILARITY_HORIZONS)
    params:
        generation_time = lambda wildcards: _generation_time_options(wildcards.analysis)
    log:
        "logs/fitness_flux/{analysis}_similarity_{horizon}d.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/forecast_similarity.py \
            --dataset {wildcards.analysis} \
            --mlr-dir mlr-estimates \
            --profile-table {input.profile} \
            --distance-matrix {input.distances} \
            --norm l1 --eps 0 \
            --slide-days {wildcards.horizon} \
            --forecast-days {wildcards.horizon} \
            {params.generation_time} \
            --detail-output {output.detail} \
            --summary-output {output.summary} 2>&1 | tee {log}
        """


rule viz_clade_distance_data:
    input:
        sequences = "fitness-flux-analysis/results/h3n2_clades_ha1_mrca.tsv",
        distances = "fitness-flux-analysis/results/h3n2_clades_distance_matrix.tsv",
        summary = "fitness-flux-analysis/results/h3n2_clades_distance_summary.json"
    output:
        data = "viz/clade-distance/data/h3n2_clades.json",
        meta = "viz/clade-distance/meta.json"
    log:
        "logs/fitness_flux/viz_clade_distance.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_clade_distance_data.py \
            --sequences {input.sequences} \
            --distances {input.distances} \
            --summary {input.summary} \
            --label H3N2 \
            --output {output.data} \
            --meta-output {output.meta} 2>&1 | tee {log}
        """


def _similarity_track_args(wildcards):
    """`--track HORIZON SUMMARY DETAIL` groups, one per forecast horizon."""
    return " ".join(
        f"--track {h} "
        f"fitness-flux-analysis/results/h3n2_clades_similarity_{h}d_summary.json "
        f"fitness-flux-analysis/results/h3n2_clades_similarity_{h}d.tsv"
        for h in SIMILARITY_HORIZONS
    )


rule viz_similarity_accuracy_data:
    input:
        summaries = expand(
            "fitness-flux-analysis/results/h3n2_clades_similarity_{horizon}d_summary.json",
            horizon=SIMILARITY_HORIZONS,
        ),
        details = expand(
            "fitness-flux-analysis/results/h3n2_clades_similarity_{horizon}d.tsv",
            horizon=SIMILARITY_HORIZONS,
        )
    output:
        data = "viz/similarity-accuracy/data/h3n2_clades.json",
        meta = "viz/similarity-accuracy/meta.json"
    params:
        tracks = _similarity_track_args
    log:
        "logs/fitness_flux/viz_similarity_accuracy.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_similarity_accuracy_data.py \
            {params.tracks} \
            --label H3N2 \
            --output {output.data} \
            --meta-output {output.meta} 2>&1 | tee {log}
        """


rule all_similarity_analysis:
    input:
        expand(
            "fitness-flux-analysis/results/{analysis}_similarity_{horizon}d_summary.json",
            analysis=SIMILARITY_ANALYSES,
            horizon=SIMILARITY_HORIZONS,
        ),
        "viz/clade-distance/data/h3n2_clades.json",
        "viz/clade-distance/meta.json",
        "viz/similarity-accuracy/data/h3n2_clades.json",
        "viz/similarity-accuracy/meta.json"
