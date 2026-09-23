#!/usr/bin/env python3
"""Build the clade-distance component's data.json.

Two panels' worth of payload, from clade_distances.py's outputs:

  * the pairwise HA1 distance matrix between clade MRCA sequences, ordered
    root-to-tip so the block structure of the tree is visible;
  * the same distances split by what the clade *labels* say about relatedness
    (ancestor/descendant, sibling, unrelated), which is the motivating panel: the
    label-based score charges the same penalty across that whole range.

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


def relationship(a, b):
    """What the clade labels alone imply about how close two clades are."""
    if b.startswith(a + ".") or a.startswith(b + "."):
        return "descendant"
    if "." in a and "." in b and a.rsplit(".", 1)[0] == b.rsplit(".", 1)[0]:
        return "sibling"
    return "unrelated"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--distances", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--label", default="H3N2")
    parser.add_argument("--output", required=True)
    parser.add_argument("--meta-output", required=True)
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
    index = {name: i for i, name in enumerate(names)}

    lookup = {}
    for row in viz_io.read_tsv(args.distances):
        lookup[(row["clade_a"], row["clade_b"])] = int(row["aa_dist"])
    matrix = [[lookup.get((a, b), None) for b in names] for a in names]

    # Upper triangle only, tagged by what the labels imply.
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pairs.append({"a": a, "b": b, "d": lookup[(a, b)], "rel": relationship(a, b)})

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
            "gene": summary["gene"],
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
                "max_distance": summary["max_distance"],
                "mean_distance": summary["mean_distance"],
                "mean_parent_child_distance": summary["mean_parent_child_distance"],
                "identical_pairs": summary["identical_pairs"],
                "aliases": summary["aliases"],
            },
        },
    )
    viz_io.write_json(
        args.meta_output,
        {"datasets": [{"id": "h3n2_clades", "label": "H3N2 clade distances"}], "default": "h3n2_clades"},
        indent=2,
    )


if __name__ == "__main__":
    main()
