#!/usr/bin/env python3
"""Per-Pango-lineage amino-acid mutation counts from the Nextclade reference tree.

Python replacement for mutation-counts-sarscov2-lineages.nb. For each tip in the
Nextclade auspice tree, walk root->leaf and count the NET amino-acid substitutions
per gene (positions whose final state differs from the Wuhan-Hu-1 reference),
rather than summing mutation events along the path.

Summing events overcounts recombinant lineages: Nextclade represents recombinants
with graft nodes (rec_parent / internal_X) carrying large overlapping/reverting
mutation lists, so the event sum balloons (e.g. XBB.1 spike 232 vs the true ~42).
The net count is the actual number of substitutions relative to Wuhan-Hu-1.

Spike (S) is split into sub-regions by amino-acid position; orf1ab = ORF1a + ORF1b;
accessory = the eight ORF/structural genes below. Counts are aggregated per Pango
lineage as the median over the tips carrying that Nextclade_pango (the pango build
has one tip per lineage, so this is normally just that tip's value).
"""

import argparse
import csv
import json
import re
import sys
from statistics import median

# Mutation string like "F486P", a deletion "H69-", or a stop "Q27*".
MUT = re.compile(r"^([A-Za-z*-])(\d+)([A-Za-z*-])$")

# Spike sub-regions by amino-acid position (matches the notebook's bounds).
SPIKE_REGIONS = {"s1": (14, 685), "s2": (686, 1273), "ntd": (14, 305), "rbd": (319, 541)}
ORF1AB_GENES = ["ORF1a", "ORF1b"]
ACCESSORY_GENES = ["ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N"]
COLUMNS = [
    "nuc_muts", "spike_muts", "s1_muts", "s2_muts",
    "ntd_muts", "rbd_muts", "orf1ab_muts", "accessory_muts",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--auspice-json", required=True, help="Nextclade auspice v2 tree (pango.json)")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def index_tree(root):
    """Flatten the tree into a list of {parent, mut, pango, tip} records."""
    nodes = []

    def walk(node, parent_idx):
        idx = len(nodes)
        attrs = node.get("node_attrs", {})
        pango = attrs.get("Nextclade_pango")
        div = attrs.get("div")
        nodes.append({
            "parent": parent_idx,
            "mut": node.get("branch_attrs", {}).get("mutations", {}),
            "pango": pango.get("value") if isinstance(pango, dict) else pango,
            "div": div.get("value") if isinstance(div, dict) else div,
            "tip": not node.get("children"),
        })
        for child in node.get("children", []):
            walk(child, idx)

    walk(root, None)
    return nodes


def net_changed_positions(nodes, leaf_idx):
    """gene -> set of positions whose final state differs from the reference, walking
    root->leaf. For each position, ref is the first mutation's 'from' and cur the last
    'to'; reversions (cur == ref) collapse away."""
    chain = []
    i = leaf_idx
    while i is not None:
        chain.append(i)
        i = nodes[i]["parent"]
    chain.reverse()

    state = {}  # gene -> {pos: [ref, cur]}
    for idx in chain:
        for gene, muts in nodes[idx]["mut"].items():
            positions = state.setdefault(gene, {})
            for mutation in muts:
                match = MUT.match(mutation)
                if not match:
                    continue
                frm, pos, to = match.group(1), int(match.group(2)), match.group(3)
                if pos in positions:
                    positions[pos][1] = to
                else:
                    positions[pos] = [frm, to]
    return {
        gene: {pos for pos, (ref, cur) in positions.items() if ref != cur}
        for gene, positions in state.items()
    }


def tip_counts(net, div):
    spike = net.get("S", set())
    counts = {
        # Nucleotide count is the tree's divergence (net substitutions from the
        # root), which excludes the masked/ambiguous sites that inflate a raw
        # per-position walk; amino-acid genes use the net per-position state.
        "nuc_muts": div if div is not None else len(net.get("nuc", set())),
        "spike_muts": len(spike),
        "orf1ab_muts": sum(len(net.get(g, set())) for g in ORF1AB_GENES),
        "accessory_muts": sum(len(net.get(g, set())) for g in ACCESSORY_GENES),
    }
    for region, (lo, hi) in SPIKE_REGIONS.items():
        counts[region + "_muts"] = sum(1 for pos in spike if lo <= pos <= hi)
    return counts


def main():
    sys.setrecursionlimit(100000)
    args = parse_args()
    with open(args.auspice_json) as handle:
        tree = json.load(handle)["tree"]
    nodes = index_tree(tree)

    per_lineage = {}
    for idx, node in enumerate(nodes):
        if not node["tip"] or not node["pango"]:
            continue
        counts = tip_counts(net_changed_positions(nodes, idx), node["div"])
        per_lineage.setdefault(node["pango"], []).append(counts)

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["lineage"] + COLUMNS)
        for lineage in sorted(per_lineage):
            tips = per_lineage[lineage]
            row = [lineage] + [f"{median([t[col] for t in tips]):g}" for col in COLUMNS]
            writer.writerow(row)

    print(f"Wrote {len(per_lineage)} lineages to {args.output}")


if __name__ == "__main__":
    main()
