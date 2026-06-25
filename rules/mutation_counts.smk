"""
Fetch the Nextclade SARS-CoV-2 reference tree and count per-Pango-lineage
amino-acid mutations (net substitutions from Wuhan-Hu-1, split by gene and by
spike sub-region). This replaces the manual mutation-counts-sarscov2-lineages.nb
step; the resulting table feeds the lineage-deltas analysis.

Counting net per-position state (rather than summing mutation events along the
root->leaf path) avoids the ~5x inflation Nextclade's recombinant graft nodes
otherwise produce for the XBB radiation. Re-run with `--forcerun fetch_nextclade_pango`
to pick up newly defined Pango lineages. Requires the `nextstrain` CLI on PATH.
"""


rule fetch_nextclade_pango:
    output:
        "data/nextclade_sarscov2_pango.json"
    params:
        url = config["mutation_counts"]["nextclade_url"],
        prefix = "data/nextclade_sarscov2_pango"
    log:
        "logs/mutation_counts/fetch.txt"
    shell:
        """
        nextstrain remote download {params.url:q} {params.prefix} 2>&1 | tee {log}
        """


rule lineage_mut_counts:
    input:
        "data/nextclade_sarscov2_pango.json"
    output:
        "mutation-counts/results/sarscov2_lineages_mut_counts.tsv"
    log:
        "logs/mutation_counts/counts.txt"
    shell:
        """
        python -u mutation-counts/count_lineage_mutations.py \
            --auspice-json {input} \
            --output {output} 2>&1 | tee {log}
        """
