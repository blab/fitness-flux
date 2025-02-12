import os, sys
import json
import Bio.Phylo
import argparse

def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.
    Assigns links back to parent nodes for the root of the tree.
    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute.
    if hasattr(node, "attr"):
        node.numdate = node.attr.get("num_date")
        node.branch_length = node.attr.get("div")

        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):
        node.branch_length = node.node_attrs.get("div")

    if root:
        node = annotate_parents_for_tree(node)

    return node

def annotate_parents_for_tree(tree):
    """Annotate each node in the given tree with its parent.
    """
    tree.root.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node

    # Return the tree.
    return tree

def collect_args():
    parser = argparse.ArgumentParser(description = "Prepare fauna FASTA for analysis")
    parser.add_argument('--json', required=True, type=str, help="JSON input file")
    parser.add_argument('--output', type=str, default=None, help="output file")

    return parser.parse_args()

if __name__=="__main__":
    params = collect_args()

    json_fh = open(params.json, "r")
    json_dict = json.load(json_fh)
    tree = json_to_tree(json_dict)

    data = []
    for n in tree.find_clades(order="postorder"):
        node_elements = {}
        node_elements["name"] = n.name
        if n.parent:
            node_elements["parent"] = n.parent.name
        else:
            node_elements["parent"] = None
        if hasattr(n, 'branch_attrs'):
            if 'mutations' in n.branch_attrs:
                node_elements["mutations"] = n.branch_attrs["mutations"]
        if hasattr(n, 'node_attrs'):
            if 'num_date' in n.node_attrs:
                if 'value' in n.node_attrs["num_date"]:
                    node_elements["num_date"] = n.node_attrs["num_date"]["value"]
            if 'clade_membership' in n.node_attrs:
                if 'value' in n.node_attrs["clade_membership"]:
                    node_elements["clade_membership"] = n.node_attrs["clade_membership"]["value"]
            if 'region' in n.node_attrs:
                if 'value' in n.node_attrs["region"]:
                    node_elements["region"] = n.node_attrs["region"]["value"]
            if 'S1_mutations' in n.node_attrs:
                if 'value' in n.node_attrs["S1_mutations"]:
                    node_elements["S1_mutations"] = n.node_attrs["S1_mutations"]["value"]
            if 'subclade' in n.node_attrs:
                if 'value' in n.node_attrs["subclade"]:
                    node_elements["subclade"] = n.node_attrs["subclade"]["value"]
            if 'Nextclade_pango' in n.node_attrs:
                if 'value' in n.node_attrs["Nextclade_pango"]:
                    node_elements["Nextclade_pango"] = n.node_attrs["Nextclade_pango"]["value"]                    
            if 'div' in n.node_attrs:
                node_elements["div"] = n.node_attrs["div"]
        print(node_elements)
        data.append(node_elements)

    with open(params.output, 'w', encoding='utf-8') as handle:
        json.dump(data, handle, indent=1, sort_keys=True)
