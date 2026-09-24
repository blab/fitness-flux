#!/usr/bin/env python3
"""Global clade x gene similarity map from an auspice tree.

Builds, once for the whole dataset, the objects the similarity-aware forecast
score needs:

  1. the reconstructed amino-acid sequence at each clade's MRCA node,
  2. the clade x (position, amino acid) indicator table ("profile table"),
  3. the pairwise clade distance matrix, in amino acids.

Originally written for the H3N2 HA1 map (--gene HA1, dot-nested subclade labels,
whole-gene distance). It also serves SARS-CoV-2, where the differences are all
parameterized rather than special-cased:

  * --clade-label first_word normalizes the tree's ``clade_membership`` values
    ("23A (XBB.1.5)") to the Nextstrain clade the MLR names ("23A");
  * --min-pos/--max-pos restrict scoring to a subunit (spike S1, aa 14-685),
    while the full gene is still reconstructed so mutations apply correctly;
  * --alias maps a merged/absent label to another clade or to the special value
    ``__root__`` (the 19A/19B "WT" merge -> Wuhan-Hu-1 root sequence);
  * a ``parent_clade`` column records each clade's nearest ancestral clade on the
    tree, so relatedness can be read from the topology rather than from label
    strings (SARS-CoV-2 clade labels are not hierarchically nested).

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
Each globally variable position (within the scored region) contributes one column
per amino acid observed at it, across all clades. A clade's row has 1 in the
column matching the amino acid it carries and 0 elsewhere. Two clades differing at
one position therefore differ in exactly two columns.

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

# Sentinel alias target: assign the root (Wuhan-Hu-1) sequence itself. Used for
# the SARS-CoV-2 "WT" bucket, which merges 19A/19B and has no MRCA node of its own.
ROOT_ALIAS = "__root__"


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
        "--gene-label",
        default=None,
        help="display label for the scored region (default: --gene), e.g. 'spike S1'",
    )
    parser.add_argument(
        "--clade-key",
        default="subclade",
        help="node_attrs key holding the clade label",
    )
    parser.add_argument(
        "--clade-label",
        choices=["exact", "first_word"],
        default="exact",
        help="'first_word' takes the label's first token, so a tree's "
        "'23A (XBB.1.5)' matches the MLR's '23A'",
    )
    parser.add_argument(
        "--min-pos",
        type=int,
        default=1,
        help="first 1-based amino-acid position scored (default 1); positions "
        "outside [min-pos, max-pos] still get their mutations applied but do not "
        "contribute to the profile table or distances (spike S1: 14)",
    )
    parser.add_argument(
        "--max-pos",
        type=int,
        default=None,
        help="last 1-based amino-acid position scored (default: gene length; spike S1: 685)",
    )
    parser.add_argument(
        "--alias",
        action="append",
        default=None,
        metavar="CLADE=TARGET",
        help="map a required clade with no MRCA node onto another clade's sequence, "
        f"or onto {ROOT_ALIAS!r} for the root sequence; repeatable",
    )
    parser.add_argument(
        "--mlr-dir",
        default=None,
        help="if set with --dataset, require every clade the MLR windows actually "
        "name; any the tree does not annotate is resolved via --alias, then to its "
        "nearest dot-delimited ancestor",
    )
    parser.add_argument("--dataset", default=None, help="e.g. h3n2_clades")
    parser.add_argument("--sequences-output", required=True)
    parser.add_argument("--distances-output", required=True)
    parser.add_argument("--profile-output", required=True)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def clade_of(node, clade_key, label_mode="exact"):
    value = node.get("node_attrs", {}).get(clade_key)
    if isinstance(value, dict):
        value = value.get("value")
    if label_mode == "first_word" and isinstance(value, str) and value:
        value = value.split()[0]
    return value if value not in EXCLUDED_CLADES else None


def mrca_sequences(tree, root_seq, gene, clade_key, label_mode="exact"):
    """{clade: sequence at its MRCA node}, tree order, and each clade's parent clade.

    A clade's MRCA is the shallowest node carrying that clade label: we walk the
    tree breadth-first from the root accumulating this gene's mutations, and the
    first time a label is seen we freeze the sequence there. Breadth-first (not
    depth-first) matters — it guarantees the shallowest node wins even if a label
    is not perfectly monophyletic.

    Alongside the sequence we record ``parents[clade]`` = the nearest ancestral
    clade label on the path from the root (or None for a basal clade). That is the
    clade's parent in the tree, which lets relatedness be read from topology rather
    than from the label string — SARS-CoV-2 Nextstrain labels are not nested.

    Mutations are accumulated into a {position: amino acid} dict rather than a
    list, so that a site mutated more than once along a path keeps only its final
    state. Deletions ('-'), stops ('*') and out-of-range positions are skipped.
    """
    length = len(root_seq)
    seqs = {}
    order = {}
    parents = {}
    # queue carries (node, inherited mutations, nearest ancestral clade label)
    queue = deque([(tree, {}, None)])
    counter = 0
    while queue:
        node, inherited, ancestor_clade = queue.popleft()
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

        clade = clade_of(node, clade_key, label_mode)
        if clade is not None and clade not in seqs:
            # apply_muts_to_root (scripts/alignment.py) expects the "A123B" form
            # and applies them in order; our dict has already resolved repeats.
            applied = [f"X{position + 1}{aa}" for position, aa in sorted(muts.items())]
            seqs[clade] = str(apply_muts_to_root(root_seq, applied))
            order[clade] = counter
            parents[clade] = ancestor_clade
            counter += 1

        # Descendants inherit this node's clade as their nearest ancestral clade
        # whenever it carries one (first occurrence or not).
        child_ancestor = clade if clade is not None else ancestor_clade
        for child in node.get("children", []) or []:
            queue.append((child, muts, child_ancestor))
    return seqs, order, parents


def parse_aliases(alias_args):
    """[`CLADE=TARGET`, ...] -> {clade: target}. TARGET may be ``__root__``."""
    aliases = {}
    for item in alias_args or []:
        if "=" not in item:
            sys.exit(f"--alias expects CLADE=TARGET, got {item!r}")
        clade, target = item.split("=", 1)
        aliases[clade.strip()] = target.strip()
    return aliases


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


def build_profile_table(seqs, positions):
    """(clades, columns, rows) one-hot over every variable scored position.

    ``positions`` is the iterable of 0-based positions eligible to be scored (the
    subunit restriction). columns are (position, amino acid) pairs; rows[i][j] is
    1 when clade i carries column j's amino acid. Two clades differing at one
    position differ in exactly two columns, so the scorer rescales by 1/2 (L1) or
    1/sqrt(2) (L2) to read the result in amino acids.
    """
    clades = sorted(seqs)

    columns = []
    for position in positions:
        alleles = sorted({seqs[clade][position] for clade in clades})
        if len(alleles) > 1:
            columns.extend((position, allele) for allele in alleles)

    rows = [
        [1 if seqs[clade][position] == allele else 0 for position, allele in columns]
        for clade in clades
    ]
    return clades, columns, rows


def hamming(a, b, positions):
    """Amino-acid differences between two sequences over the scored positions."""
    return sum(1 for p in positions if a[p] != b[p])


def main():
    args = parse_args()
    gene_label = args.gene_label or args.gene
    explicit_aliases = parse_aliases(args.alias)

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
    length = len(root_seq)

    # Scored region (subunit restriction). Positions outside it still had their
    # mutations applied during reconstruction, but do not contribute to columns or
    # distances. Reported in 1-based coordinates; defaults span the whole gene.
    min_pos = max(1, args.min_pos)
    max_pos = min(length, args.max_pos if args.max_pos is not None else length)
    positions = list(range(min_pos - 1, max_pos))  # 0-based scored positions
    restricted = (min_pos > 1) or (max_pos < length)

    seqs, order, parents = mrca_sequences(
        tree_json["tree"], root_seq, args.gene, args.clade_key, args.clade_label
    )
    ff_io.log(
        f"Reconstructed {args.gene} at the MRCA of {len(seqs)} clades ({length} aa); "
        f"scoring {gene_label} positions {min_pos}-{max_pos} ({len(positions)} aa)"
    )

    # Resolve any required clade the tree never annotates: explicit --alias first
    # (including the __root__ sentinel), then the dot-nomenclature ancestor fallback.
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
            target = explicit_aliases.get(clade)
            if target == ROOT_ALIAS:
                seqs[clade] = str(root_seq)
                order[clade] = -1  # basal: sort ahead of every reconstructed clade
                parents[clade] = None
                aliases[clade] = ROOT_ALIAS
                ff_io.log(f"  {clade}: aliased to the root sequence")
                continue
            if target and target in seqs:
                aliases[clade] = target
                ff_io.log(f"  {clade}: no MRCA node, using --alias target {target}")
                continue
            fallback = target or nearest_ancestor(clade, seqs)
            if fallback and fallback in seqs:
                aliases[clade] = fallback
                ff_io.log(f"  {clade}: no MRCA node in tree, using ancestor {fallback}")
            else:
                missing.append(clade)
                ff_io.log(f"  {clade}: no MRCA node and no alias/ancestor in tree - EXCLUDED")

    # Aliased clades become first-class rows carrying their target's sequence, so
    # every downstream consumer covers them without special-casing. They land at
    # distance 0 from that target, which is exactly what the scorer's identity
    # columns exist to separate.
    for alias, target in aliases.items():
        if target == ROOT_ALIAS:
            continue  # already assigned the root sequence above
        seqs[alias] = seqs[target]
        order[alias] = order[target]
        parents[alias] = target  # its target is its nearest present relative

    clades, columns, rows = build_profile_table(seqs, positions)
    variable_positions = sorted({position for position, _ in columns})
    ff_io.log(
        f"{len(columns)} profile columns over {len(variable_positions)} variable "
        f"{gene_label} positions"
    )

    with open(args.sequences_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["clade", "tree_order", "parent_clade", "n_subs_from_root", args.gene.lower() + "_seq"]
        )
        for clade in clades:
            writer.writerow(
                [
                    clade,
                    order[clade],
                    parents.get(clade) or "",
                    hamming(seqs[clade], root_seq, positions),
                    seqs[clade],
                ]
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
                d = hamming(seqs[a], seqs[b], positions)
                distances[(a, b)] = d
                writer.writerow([a, b, d])
    ff_io.log(
        f"Wrote {len(clades)}x{len(clades)} distance matrix to {args.distances_output}"
    )

    offdiag = [distances[(a, b)] for a in clades for b in clades if a < b]
    offdiag.sort()
    # Parent-child diagnostic (summary only). Prefer the dot-nomenclature test so
    # the published H3N2 number is unchanged; when labels are not nested (SARS-CoV-2)
    # no dot pair matches, so fall back to the tree topology in parent_clade. The
    # figure's relatedness panel always uses the topology (parent_clade).
    parent_child = [
        (a, b, distances[(a, b)])
        for a in clades
        for b in clades
        if b.startswith(a + ".") and b.count(".") == a.count(".") + 1
    ]
    if not parent_child:
        clade_set = set(clades)
        parent_child = [
            (parents[b], b, distances[(parents[b], b)])
            for b in clades
            if parents.get(b) in clade_set
        ]
    identical = [(a, b) for a in clades for b in clades if a < b and distances[(a, b)] == 0]

    summary = {
        "gene": args.gene,
        "gene_label": gene_label,
        "clade_key": args.clade_key,
        "n_clades": len(clades),
        "gene_length": length,
        "scored_region": [min_pos, max_pos] if restricted else None,
        "n_scored_positions": len(positions),
        "n_variable_positions": len(variable_positions),
        "n_profile_columns": len(columns),
        "variable_positions": [position + 1 for position in variable_positions],
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
