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

import json


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


rule merge_clades:
    """Relabel (merge) clade values in the provisioned metadata per the per-virus
    `merge_clades` config (e.g. SARS-CoV-2 19A + 19B -> WT), writing the file each
    sarscov2 dataset's local_metadata points at. Viruses with no merge_clades pass
    through unchanged. Reads the existing subset file, so no re-download."""
    input:
        "data/{virus}_subset_metadata.tsv.zst"
    output:
        "data/{virus}_merged_metadata.tsv.zst"
    wildcard_constraints:
        virus = "sarscov2|h3n2"
    params:
        merges = lambda w: json.dumps(config["provision"][w.virus].get("merge_clades", {}))
    log:
        "logs/provision/{virus}_merge.txt"
    shell:
        """
        set -euo pipefail
        ( zstd -dc {input} \
            | python scripts/merge-clades.py --merges {params.merges:q} \
            | zstd -c > {output} ) 2> {log}
        """


rule all_provision_metadata:
    input:
        sorted({config[dataset]["local_metadata"] for dataset in config["datasets"]})
