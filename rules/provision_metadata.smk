"""
Provision the per-virus input metadata that the rest of the workflow consumes.

Each dataset block's `local_metadata` points at data/{virus}_subset_metadata.tsv.zst.
This rule streams the full Nextstrain metadata from S3, decompresses it, subsets
to the columns the pipeline uses, and recompresses to that file — replacing the
manual aws/zstd/tsv-select steps that used to live in the README. The source URL,
decompressor, and column set are per-virus config under `provision`.

The rule has no `input:`, so Snakemake provisions each file once and never
re-downloads it on its own; refresh with `--forcerun provision_metadata` or by
deleting the file. (data/ is gitignored and untouched by the clean rules.)

Requires `aws`, `zstd`, `xz`, and `tsv-select` on PATH; H3N2 reads a private
bucket and needs AWS credentials.
"""

rule provision_metadata:
    output:
        "data/{virus}_subset_metadata.tsv.zst"
    wildcard_constraints:
        virus = "sarscov2|h3n2"
    params:
        url = lambda w: config["provision"][w.virus]["metadata_url"],
        decompress = lambda w: config["provision"][w.virus]["decompress"],
        columns = lambda w: config["provision"][w.virus]["columns"]
    log:
        "logs/provision/{virus}.txt"
    shell:
        """
        set -euo pipefail
        ( aws s3 cp {params.url} - \
            | {params.decompress} \
            | tsv-select -H -f {params.columns} \
            | zstd -c > {output} ) 2> {log}
        """

rule all_provision_metadata:
    input:
        expand("data/{virus}_subset_metadata.tsv.zst", virus=["sarscov2", "h3n2"])
