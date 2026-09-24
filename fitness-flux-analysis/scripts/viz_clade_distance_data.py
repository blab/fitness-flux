#!/usr/bin/env python3
"""Build the clade-distance component's data.json.

Two panels' worth of payload, from clade_distances.py's outputs:

  * the pairwise distance matrix between clade MRCA sequences (HA1 for H3N2,
    spike S1 for SARS-CoV-2), ordered root-to-tip so the block structure of the
    tree is visible;
  * the same distances split by how related two clades are (ancestor/descendant,
    sibling, unrelated), which is the motivating panel: the label-based score
    charges the same penalty across that whole range.

Relatedness is determined by --relatedness:

  * ``label`` (default): from the dot-nested clade label strings (H3N2 subclades,
    e.g. J.2 is an ancestor of J.2.4);
  * ``tree``: from the tree topology in the sequences file's ``parent_clade``
    column (SARS-CoV-2 Nextstrain clades like 21K/22B are not nested, so the
    label string carries no ancestry — the tree does).

The component renders these; it does no joining or classification of its own.
"""

import argparse
import json

import viz_io

# Sequential single hue, light to dark, for the distance heatmap (dataviz skill:
# magnitude gets one hue, never a rainbow). Anchored on the same blue as the MLR
# series so the figure set reads as one system.
SEQUENTIAL = ["#eef4fc", "#cfe0f6", "#a8c8ee", "#7aa9e2", "#4a86d4", "#2a78d6", "#1a4f93"]

RELATIONSHIPS = [
    {"key": "descendant", "label": "Ancestor/descendant", "color": "#2a78d6"},
    {"key": "sibling", "label": "Sibling", "color": "#7aa9e2"},
    {"key": "unrelated", "label": "Unrelated label", "color": "#b8b6ae"},
]


def label_relationship(a, b):
    """What the dot-nested clade labels alone imply about how close two clades are."""
    if b.startswith(a + ".") or a.startswith(b + "."):
        return "descendant"
    if "." in a and "." in b and a.rsplit(".", 1)[0] == b.rsplit(".", 1)[0]:
        return "sibling"
    return "unrelated"


def ancestors_of(clade, parent):
    """The set of ancestral clades on `clade`'s path to the root, via parent_clade."""
    seen = set()
    cur = parent.get(clade)
    while cur and cur not in seen:
        seen.add(cur)
        cur = parent.get(cur)
    return seen


def tree_relationship(a, b, parent):
    """Relatedness read from the tree topology (the parent_clade column)."""
    if a in ancestors_of(b, parent) or b in ancestors_of(a, parent):
        return "descendant"
    pa, pb = parent.get(a), parent.get(b)
    if pa is not None and pa == pb:
        return "sibling"
    return "unrelated"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--label", default="H3N2")
    parser.add_argument(
        "--relatedness",
        choices=["label", "tree"],
        default="label",
        help="classify clade pairs by dot-nested label strings (H3N2) or by the "
        "tree topology in parent_clade (SARS-CoV-2)",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--meta-output",
        default=None,
        help="optional; the similarity components' meta.json is normally written "
        "once by the viz_similarity_meta rule",
    )
    args = parser.parse_args()

    with open(args.summary) as handle:
        summary = json.load(handle)

    # Order clades root-to-tip by divergence from the root, so the heatmap's
    # block structure follows the tree rather than the alphabet.
    records = viz_io.read_tsv(args.sequences)
    clades = sorted(
        ({"clade": r["clade"], "n_subs": int(r["n_subs_from_root"]), "order": int(r["tree_order"])}
         for r in records),
        key=lambda r: (r["n_subs"], r["order"], r["clade"]),
    )
    names = [r["clade"] for r in clades]
    parent = {r["clade"]: (r.get("parent_clade") or None) for r in records}

    def classify(a, b):
        return tree_relationship(a, b, parent) if args.relatedness == "tree" else label_relationship(a, b)

    lookup = {}
    for row in viz_io.read_tsv(args.distances):
        lookup[(row["clade_a"], row["clade_b"])] = int(row["aa_dist"])
    matrix = [[lookup.get((a, b), None) for b in names] for a in names]

    # Upper triangle only, tagged by how related the two clades are.
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pairs.append({"a": a, "b": b, "d": lookup[(a, b)], "rel": classify(a, b)})

    by_rel = {}
    for entry in pairs:
        by_rel.setdefault(entry["rel"], []).append(entry["d"])
    stats = {
        key: {
            "n": len(values),
            "mean": round(sum(values) / len(values), 2),
            "min": min(values),
            "max": max(values),
        }
        for key, values in by_rel.items()
    }

    viz_io.write_json(
        args.output,
        {
            "label": args.label,
            "gene": summary.get("gene_label", summary["gene"]),
            "clades": clades,
            "names": names,
            "matrix": matrix,
            "pairs": pairs,
            "relationships": RELATIONSHIPS,
            "stats": stats,
            "sequential": SEQUENTIAL,
            "summary": {
                "n_clades": summary["n_clades"],
                "n_variable_positions": summary["n_variable_positions"],
                "n_profile_columns": summary["n_profile_columns"],
                "gene_length": summary["gene_length"],
                "scored_region": summary.get("scored_region"),
                "max_distance": summary["max_distance"],
                "mean_distance": summary["mean_distance"],
                "mean_parent_child_distance": summary["mean_parent_child_distance"],
                "identical_pairs": summary["identical_pairs"],
                "aliases": summary["aliases"],
            },
        },
    )
    if args.meta_output:
        viz_io.write_json(
            args.meta_output,
            {"datasets": [{"id": "h3n2_clades", "label": "H3N2 clade distances"}], "default": "h3n2_clades"},
            indent=2,
        )


if __name__ == "__main__":
    main()
