#!/usr/bin/env python3
"""Global clade x HA1 similarity map from a seasonal-flu auspice tree.

Builds, once for the whole dataset, the objects the similarity-aware forecast
score needs:

  1. the reconstructed HA1 amino-acid sequence at each clade's MRCA node,
  2. the clade x (position, amino acid) indicator table ("profile table"),
  3. the pairwise clade distance matrix, in amino acids.

Why the MRCA node and not a consensus
-------------------------------------
A clade's consensus sequence is a mixture over whatever descendants happen to be
assigned to it, so it shifts whenever the upstream collapse threshold changes.
The MRCA node sequence is a property of the tree alone and is stable.

Why a single GLOBAL map
-----------------------
The forecast scorer compares a forecast made over window W's clade set against
truth observed over window W+1's clade set. Those sets differ. Rather than
remapping one onto the other (which is what drops or mislabels newly emerged
clades in the first place), both are embedded in one shared space spanning every
clade in the tree. Nothing has to be remapped, renormalized onto a common label
set, or dropped.

The encoding
------------
Each globally variable HA1 position contributes one column per amino acid
observed at it, across all clades. A clade's row has 1 in the column matching
the amino acid it carries and 0 elsewhere. Two clades differing at one position
therefore differ in exactly two columns.

The table is written as plain 0/1; the scorer applies the scale its norm needs
(1/2 under L1, 1/sqrt(2) under L2) so that the distance between two pure clades
comes out in whole amino acids either way.

(One column per allele, rather than one column per non-reference allele, is what
makes this consistent. Reference-allele encoding gives a biallelic site one
column and a triallelic site two, so a single amino-acid difference would count
differently depending on how many other clades happen to vary at that site.)

Identity columns
----------------
A further column per clade, carrying `eps` on the diagonal, is appended by the
scorer (not here). Those guarantee distinct clades never collapse onto the same
row. They are not cosmetic: in the H3N2 tree at least one clade pair has
identical reconstructed HA1 MRCA sequences.
"""

import argparse
import csv
import json
import os
import sys
from collections import deque

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
from alignment import apply_muts_to_root, load_json  # noqa: E402

import ff_io  # noqa: E402

# Clade labels that carry no reconstructable sequence: "unassigned" is the
# tree's catch-all for tips outside any named clade, "other" is the analysis
# pipeline's collapsed rare-clade bucket.
EXCLUDED_CLADES = {"unassigned", "other", "", None}


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="auspice tree JSON (path or URL)")
    parser.add_argument(
        "--root",
        default=None,
        help="root-sequence JSON; omit when the tree JSON embeds 'root_sequence'",
    )
    parser.add_argument("--gene", default="HA1", help="gene to reconstruct")
    parser.add_argument(
        "--clade-key",
        default="subclade",
        help="node_attrs key holding the clade label",
    )
    parser.add_argument(
        "--mlr-dir",
        default=None,
        help="if set with --dataset, require every clade the MLR windows actually "
        "name; any the tree does not annotate is resolved to its nearest ancestor",
    )
    parser.add_argument("--dataset", default=None, help="e.g. h3n2_clades")
    parser.add_argument("--sequences-output", required=True)
    parser.add_argument("--distances-output", required=True)
    parser.add_argument("--profile-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def clade_of(node, clade_key):
    value = node.get("node_attrs", {}).get(clade_key)
    if isinstance(value, dict):
        value = value.get("value")
    return value if value not in EXCLUDED_CLADES else None


def mrca_sequences(tree, root_seq, gene, clade_key):
    """{clade: HA1 sequence at its MRCA node}, plus the root-to-node substitution count.

    A clade's MRCA is the shallowest node carrying that clade label: we walk the
    tree breadth-first from the root accumulating this gene's mutations, and the
    first time a label is seen we freeze the sequence there. Breadth-first (not
    depth-first) matters — it guarantees the shallowest node wins even if a label
    is not perfectly monophyletic.

    Mutations are accumulated into a {position: amino acid} dict rather than a
    list, so that a site mutated more than once along a path keeps only its final
    state. Deletions ('-'), stops ('*') and out-of-range positions are skipped.
    """
    length = len(root_seq)
    seqs = {}
    order = {}
    queue = deque([(tree, {})])
    counter = 0
    while queue:
        node, inherited = queue.popleft()
        muts = dict(inherited)
        for mut in (node.get("branch_attrs", {}).get("mutations", {}) or {}).get(gene, []) or []:
            if "*" in mut or "-" in mut:
                continue
            try:
                position = int(mut[1:-1]) - 1
            except ValueError:
                continue
            if 0 <= position < length:
                muts[position] = mut[-1]

        clade = clade_of(node, clade_key)
        if clade is not None and clade not in seqs:
            # apply_muts_to_root (scripts/alignment.py) expects the "A123B" form
            # and applies them in order; our dict has already resolved repeats.
            applied = [f"X{position + 1}{aa}" for position, aa in sorted(muts.items())]
            seqs[clade] = str(apply_muts_to_root(root_seq, applied))
            order[clade] = counter
            counter += 1

        for child in node.get("children", []) or []:
            queue.append((child, muts))
    return seqs, order


def nearest_ancestor(clade, available):
    """Nearest dot-delimited ancestor of `clade` present in `available`, else None.

    Used only as a fallback for a clade the tree never annotates as its own node
    (H3N2's B.1.2 is one: the label exists in the metadata but no node carries
    it). Flu subclade nomenclature is strictly nested by dots, so stripping the
    last component walks up the hierarchy.
    """
    parts = clade.split(".")
    for cut in range(len(parts) - 1, 0, -1):
        candidate = ".".join(parts[:cut])
        if candidate in available:
            return candidate
    return None


def build_profile_table(seqs):
    """(clades, columns, rows) one-hot over every globally variable position.

    columns are (position, amino acid) pairs; rows[i][j] is 1 when clade i
    carries column j's amino acid. Two clades differing at one position differ in
    exactly two columns, so the scorer rescales by 1/2 (L1) or 1/sqrt(2) (L2) to
    read the result in amino acids.
    """
    clades = sorted(seqs)
    length = len(next(iter(seqs.values())))

    columns = []
    for position in range(length):
        alleles = sorted({seqs[clade][position] for clade in clades})
        if len(alleles) > 1:
            columns.extend((position, allele) for allele in alleles)

    rows = [
        [1 if seqs[clade][position] == allele else 0 for position, allele in columns]
        for clade in clades
    ]
    return clades, columns, rows


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def main():
    args = parse_args()

    tree_json = load_json(args.tree)
    if args.root:
        root_json = load_json(args.root)
    elif "root_sequence" in tree_json:
        root_json = tree_json["root_sequence"]
    else:
        sys.exit("No --root given and the tree JSON has no embedded 'root_sequence'")

    if args.gene not in root_json:
        sys.exit(f"Gene {args.gene!r} not in root sequence (have: {sorted(root_json)})")
    root_seq = root_json[args.gene]

    seqs, order = mrca_sequences(tree_json["tree"], root_seq, args.gene, args.clade_key)
    ff_io.log(f"Reconstructed {args.gene} at the MRCA of {len(seqs)} clades ({len(root_seq)} aa)")

    # Resolve any required clade the tree never annotates onto its nearest ancestor.
    aliases = {}
    missing = []
    if args.mlr_dir and args.dataset:
        required = set()
        for timepoint in ff_io.seasonal_timepoints(args.mlr_dir, args.dataset):
            mlr = ff_io.load_mlr(args.mlr_dir, args.dataset, timepoint)
            required.update(mlr["metadata"]["variants"])
        ff_io.log(f"{len(required)} distinct clades named across {args.dataset} windows")
        for clade in sorted(required):
            if clade in seqs or clade in EXCLUDED_CLADES:
                continue
            fallback = nearest_ancestor(clade, seqs)
            if fallback:
                aliases[clade] = fallback
                ff_io.log(f"  {clade}: no MRCA node in tree, using ancestor {fallback}")
            else:
                missing.append(clade)
                ff_io.log(f"  {clade}: no MRCA node and no ancestor in tree - EXCLUDED")

    # Aliased clades become first-class rows carrying their ancestor's sequence,
    # so every downstream consumer covers them without special-casing. They land
    # at distance 0 from that ancestor, which is exactly what the scorer's
    # identity columns exist to separate.
    for alias, target in aliases.items():
        seqs[alias] = seqs[target]
        order[alias] = order[target]

    clades, columns, rows = build_profile_table(seqs)
    positions = sorted({position for position, _ in columns})
    ff_io.log(f"{len(columns)} profile columns over {len(positions)} variable {args.gene} positions")

    with open(args.sequences_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["clade", "tree_order", "n_subs_from_root", args.gene.lower() + "_seq"])
        for clade in clades:
            writer.writerow(
                [clade, order[clade], hamming(seqs[clade], root_seq), seqs[clade]]
            )
    ff_io.log(f"Wrote {len(clades)} MRCA sequences to {args.sequences_output}")

    with open(args.profile_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["clade"] + [f"{position + 1}{allele}" for position, allele in columns])
        for clade, row in zip(clades, rows):
            writer.writerow([clade] + row)
    ff_io.log(f"Wrote profile table to {args.profile_output}")

    # Full matrix including the diagonal and both orientations, so consumers can
    # index it either way without worrying about ordering.
    distances = {}
    with open(args.distances_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["clade_a", "clade_b", "aa_dist"])
        for a in clades:
            for b in clades:
                d = hamming(seqs[a], seqs[b])
                distances[(a, b)] = d
                writer.writerow([a, b, d])
    ff_io.log(
        f"Wrote {len(clades)}x{len(clades)} distance matrix to {args.distances_output}"
    )

    offdiag = [distances[(a, b)] for a in clades for b in clades if a < b]
    offdiag.sort()
    parent_child = [
        (a, b, distances[(a, b)])
        for a in clades
        for b in clades
        if b.startswith(a + ".") and b.count(".") == a.count(".") + 1
    ]
    identical = [(a, b) for a in clades for b in clades if a < b and distances[(a, b)] == 0]

    summary = {
        "gene": args.gene,
        "clade_key": args.clade_key,
        "n_clades": len(clades),
        "gene_length": len(root_seq),
        "n_variable_positions": len(positions),
        "n_profile_columns": len(columns),
        "variable_positions": [position + 1 for position in positions],
        "aliases": aliases,
        "missing": missing,
        "max_distance": offdiag[-1] if offdiag else 0,
        "median_distance": offdiag[len(offdiag) // 2] if offdiag else 0,
        "mean_distance": round(sum(offdiag) / len(offdiag), 3) if offdiag else 0,
        "mean_parent_child_distance": (
            round(sum(d for _, _, d in parent_child) / len(parent_child), 3)
            if parent_child else None
        ),
        "n_parent_child_pairs": len(parent_child),
        "identical_pairs": [list(pair) for pair in identical],
    }
    with open(args.summary_output, "w") as handle:
        json.dump(summary, handle, indent=2)
    ff_io.log(
        f"max {summary['max_distance']} aa, mean {summary['mean_distance']} aa, "
        f"parent-child mean {summary['mean_parent_child_distance']} aa, "
        f"{len(identical)} identical pair(s)"
    )


if __name__ == "__main__":
    main()
