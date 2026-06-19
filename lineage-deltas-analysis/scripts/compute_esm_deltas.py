#!/usr/bin/env python3
"""Compute per-branch ESM embedding distance vs fitness change.

Port of the "ESM correlation" section of ``lineage-deltas.nb``. Each lineage has
a protein-language-model embedding (ESM-2, 650M and 3B parameters, pretrained and
fine-tuned). Along a parent->child branch we take the Euclidean distance between
the endpoint embeddings and pair it with the per-season change in log fitness.

Embedding files are large and xz-compressed; only the embeddings for lineages
that actually appear on a branch are retained.
"""

import argparse
import csv
import lzma

import numpy as np

import ld_io

EMBEDDINGS = {
    "esm_650M_pretrained": "embeddings_650M_pretrained.tsv.xz",
    "esm_650M_fine_tuned": "embeddings_650M_fine_tuned.tsv.xz",
    "esm_3B_pretrained": "embeddings_3B_pretrained.tsv.xz",
    "esm_3B_fine_tuned": "embeddings_3B_fine_tuned.tsv.xz",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--mlr-dir", default="mlr-estimates")
    parser.add_argument("--seqcounts-dir", default="sequence-counts")
    parser.add_argument(
        "--source-data", default="lineage-deltas-analysis/source-data"
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def load_embeddings(path, wanted):
    """name -> np.array embedding, restricted to ``wanted`` lineages."""
    embeddings = {}
    with lzma.open(path, "rt") as handle:
        for line in handle:
            name, _, rest = line.partition("\t")
            if name in wanted:
                embeddings[name] = np.array(rest.split(), dtype=float)
    return embeddings


def main():
    args = parse_args()
    branches = ld_io.branches(args.mlr_dir, args.seqcounts_dir)
    wanted = {b["parent"] for b in branches} | {b["child"] for b in branches}

    distances = {key: {} for key in EMBEDDINGS}
    for key, filename in EMBEDDINGS.items():
        path = f"{args.source_data}/{filename}"
        embeddings = load_embeddings(path, wanted)
        for branch in branches:
            child, parent = branch["child"], branch["parent"]
            if child in embeddings and parent in embeddings:
                distances[key][(branch["timepoint"], parent, child)] = float(
                    np.linalg.norm(embeddings[child] - embeddings[parent])
                )
        ld_io.log(f"{key}: {len(distances[key])} branches with embeddings")

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        columns = list(EMBEDDINGS)
        writer.writerow(
            ["timepoint", "parent", "child", "delta_log_fitness"] + columns
        )
        for branch in branches:
            key = (branch["timepoint"], branch["parent"], branch["child"])
            values = [distances[col].get(key) for col in columns]
            if all(v is None for v in values):
                continue
            writer.writerow(
                [
                    branch["timepoint"],
                    branch["parent"],
                    branch["child"],
                    f"{branch['delta_log_fitness']:.6f}",
                ]
                + ["" if v is None else f"{v:.6f}" for v in values]
            )
    ld_io.log(f"Wrote ESM distance deltas to {args.output}")


if __name__ == "__main__":
    main()
