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

# Display labels for the similarity figures/tables, per dataset.
SIMILARITY_LABELS = {"h3n2_clades": "H3N2", "sarscov2_clades": "SARS-CoV-2"}


def _clade_tree(analysis, field):
    return config["clade_distance"]["trees"][analysis][field]


def _clade_tree_entry(analysis):
    return config["clade_distance"]["trees"][analysis]


def _clade_tree_path(analysis):
    """The tree JSON for a dataset: a local file (SARS-CoV-2, fetched by
    download_ncov_tree) when `local` is set, else the download_clade_tree output."""
    return _clade_tree_entry(analysis).get("local", f"data/{analysis}_clade_tree.json")


def _clade_distances_inputs(wildcards):
    """Tree (+ optional root sidecar) and the MLR windows the map must cover."""
    entry = _clade_tree_entry(wildcards.analysis)
    inputs = {
        "tree": _clade_tree_path(wildcards.analysis),
        "mlr": _forecast_season_inputs(wildcards),
    }
    if entry.get("root"):
        inputs["root"] = entry["root"]
    return inputs


def _clade_distance_args(wildcards):
    """Assemble clade_distances.py's per-dataset flags from the tree config.

    H3N2 gets only --gene/--clade-key (its defaults); SARS-CoV-2 adds the root
    sidecar, the display label, first-word label normalization, the S1 position
    window, and the WT->root alias."""
    entry = _clade_tree_entry(wildcards.analysis)
    parts = [f"--gene {entry['gene']}", f"--clade-key {entry['clade_key']}"]
    if entry.get("root"):
        parts.append(f"--root {entry['root']}")
    if entry.get("gene_label"):
        parts.append(f"--gene-label {entry['gene_label']!r}")
    if entry.get("clade_label"):
        parts.append(f"--clade-label {entry['clade_label']}")
    if entry.get("min_pos") is not None:
        parts.append(f"--min-pos {entry['min_pos']}")
    if entry.get("max_pos") is not None:
        parts.append(f"--max-pos {entry['max_pos']}")
    for clade, target in (entry.get("aliases") or {}).items():
        parts.append(f"--alias {clade}={target}")
    return " ".join(parts)


def _similarity_relatedness(analysis):
    """How the clade-distance figure classifies relatedness: from the tree
    topology when labels are not hierarchically nested (SARS-CoV-2, first_word),
    else from the dot-nested label strings (H3N2)."""
    return "tree" if _clade_tree_entry(analysis).get("clade_label") == "first_word" else "label"


def _similarity_meta_datasets():
    """`id:Label,...` for viz_meta.py, over the datasets that have a tree."""
    return ",".join(f"{a}:{SIMILARITY_LABELS.get(a, a)} clades" for a in SIMILARITY_ANALYSES)


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
        unpack(_clade_distances_inputs)
    output:
        sequences = "fitness-flux-analysis/results/{analysis}_ha1_mrca.tsv",
        distances = "fitness-flux-analysis/results/{analysis}_distance_matrix.tsv",
        profile = "fitness-flux-analysis/results/{analysis}_profile_table.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_distance_summary.json"
    params:
        args = _clade_distance_args
    log:
        "logs/fitness_flux/{analysis}_clade_distances.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/clade_distances.py \
            --tree {input.tree} \
            {params.args} \
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
            --endpoint 0 --endpoint 30 --endpoint 90 --endpoint 180 --endpoint 365 \
            {params.generation_time} \
            --detail-output {output.detail} \
            --summary-output {output.summary} 2>&1 | tee {log}
        """


rule viz_clade_distance_data:
    """Per-dataset clade-distance payload. Relatedness is read from the tree
    topology for non-nested labels (SARS-CoV-2) and from the label strings for
    dot-nested labels (H3N2). meta.json is written once by viz_similarity_meta."""
    input:
        sequences = "fitness-flux-analysis/results/{analysis}_ha1_mrca.tsv",
        distances = "fitness-flux-analysis/results/{analysis}_distance_matrix.tsv",
        summary = "fitness-flux-analysis/results/{analysis}_distance_summary.json"
    output:
        data = "viz/clade-distance/data/{analysis}.json"
    wildcard_constraints:
        analysis = "|".join(SIMILARITY_ANALYSES)
    params:
        label = lambda w: SIMILARITY_LABELS.get(w.analysis, w.analysis),
        relatedness = lambda w: _similarity_relatedness(w.analysis)
    log:
        "logs/fitness_flux/{analysis}_viz_clade_distance.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_clade_distance_data.py \
            --sequences {input.sequences} \
            --distances {input.distances} \
            --summary {input.summary} \
            --label {params.label:q} \
            --relatedness {params.relatedness} \
            --output {output.data} 2>&1 | tee {log}
        """


def _similarity_track_args(wildcards):
    """`--track HORIZON SUMMARY DETAIL` groups, one per forecast horizon."""
    return " ".join(
        f"--track {h} "
        f"fitness-flux-analysis/results/{wildcards.analysis}_similarity_{h}d_summary.json "
        f"fitness-flux-analysis/results/{wildcards.analysis}_similarity_{h}d.tsv"
        for h in SIMILARITY_HORIZONS
    )


rule viz_similarity_accuracy_data:
    """Per-dataset similarity-accuracy payload. meta.json is written once by
    viz_similarity_meta."""
    input:
        summaries = expand(
            "fitness-flux-analysis/results/{{analysis}}_similarity_{horizon}d_summary.json",
            horizon=SIMILARITY_HORIZONS,
        ),
        details = expand(
            "fitness-flux-analysis/results/{{analysis}}_similarity_{horizon}d.tsv",
            horizon=SIMILARITY_HORIZONS,
        )
    output:
        data = "viz/similarity-accuracy/data/{analysis}.json"
    wildcard_constraints:
        analysis = "|".join(SIMILARITY_ANALYSES)
    params:
        tracks = _similarity_track_args,
        label = lambda w: SIMILARITY_LABELS.get(w.analysis, w.analysis),
        gene_label = lambda w: _clade_tree_entry(w.analysis).get("gene_label", "HA1")
    log:
        "logs/fitness_flux/{analysis}_viz_similarity_accuracy.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_similarity_accuracy_data.py \
            {params.tracks} \
            --label {params.label:q} \
            --gene-label {params.gene_label:q} \
            --output {output.data} 2>&1 | tee {log}
        """


rule viz_similarity_meta:
    """Each similarity component's meta.json: the datasets it carries (those with a
    clade_distance tree), with SARS-CoV-2 as the default to match the other components."""
    output:
        clade_distance = "viz/clade-distance/meta.json",
        similarity_accuracy = "viz/similarity-accuracy/meta.json"
    params:
        datasets = _similarity_meta_datasets(),
        default = SIMILARITY_ANALYSES[0] if SIMILARITY_ANALYSES else ""
    log:
        "logs/fitness_flux/viz_similarity_meta.txt"
    shell:
        """
        python -u fitness-flux-analysis/scripts/viz_meta.py \
            --datasets {params.datasets:q} --default {params.default} \
            --output {output.clade_distance} 2>&1 | tee {log}
        python -u fitness-flux-analysis/scripts/viz_meta.py \
            --datasets {params.datasets:q} --default {params.default} \
            --output {output.similarity_accuracy} 2>&1 | tee -a {log}
        """


rule all_similarity_analysis:
    input:
        expand(
            "fitness-flux-analysis/results/{analysis}_similarity_{horizon}d_summary.json",
            analysis=SIMILARITY_ANALYSES,
            horizon=SIMILARITY_HORIZONS,
        ),
        expand("viz/clade-distance/data/{analysis}.json", analysis=SIMILARITY_ANALYSES),
        expand("viz/similarity-accuracy/data/{analysis}.json", analysis=SIMILARITY_ANALYSES),
        "viz/clade-distance/meta.json",
        "viz/similarity-accuracy/meta.json"
