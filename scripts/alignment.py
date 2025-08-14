# Code from Katie Kistler
# MIT license
"""
Given a tree.json and root-sequence.json file, finds the sequences of
each node in the tree and outputs a FASTA file with these node sequences.
If a gene is specified, the sequences will be the AA sequence of that gene
at that node. If 'nuc' is specified, the whole genome nucleotide sequence
at the node will be output. (this is default if no gene is specified).
The FASTA header is the node's name in the tree.json
"""
import argparse
import json
import requests
from augur.utils import json_to_tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from urllib.parse import urlparse

def apply_muts_to_root(root_seq, list_of_muts):
    """
    Apply a list of mutations to the root sequence
    to find the sequence at a given node. The list of mutations
    is ordered from root to node, so multiple mutations at the
    same site will correctly overwrite each other
    """
    # make the root sequence mutable
    root_plus_muts = MutableSeq(root_seq)

    # apply all mutations to root sequence
    for mut in list_of_muts:
        # subtract 1 to deal with biological numbering vs python
        mut_site = int(mut[1:-1])-1
        # get the nucleotide that the site was mutated to
        mutation = mut[-1]
        # apply mutation
        root_plus_muts[mut_site] = mutation

    return root_plus_muts

def load_json(source):
    """Load JSON data from a URL or local file."""
    parsed = urlparse(source)
    if parsed.scheme in ('http', 'https'):
        return requests.get(source, headers={"accept": "application/json"}).json()
    else:
        with open(source, 'r') as f:
            return json.load(f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--gene", default="nuc",
                        help="Name of gene to return AA sequences for. 'nuc' will return full genome nucleotide seq")
    parser.add_argument("--tree", required=True,
                        help="URL or local path for the tree.json file")
    parser.add_argument("--root", required=True,
                        help="URL or local path for the root-sequence.json file")
    parser.add_argument("--output", type=str, default="alignment.fasta", help="Output FASTA file for sequences")
    parser.add_argument("--tips-only", action="store_true", help="If set, only include tip sequences (leaf nodes)")

    args = parser.parse_args()

    # Load tree and root JSON from either URLs or local files
    tree_json = load_json(args.tree)
    root_json = load_json(args.root)

    # Convert tree JSON to Bio.phylo format
    tree = json_to_tree(tree_json)

    # Get the nucleotide sequence of the root
    root_seq = root_json[args.gene]

    # Initialize list to store sequence records for each node
    sequence_records = []

    # Iterate over tip nodes only if --tips-only is set, otherwise iterate over all nodes
    nodes = list(tree.get_terminals() if args.tips_only else tree.find_clades())

    for node in tqdm(nodes, desc="Processing nodes", unit="node"):

        # Get path back to the root
        path = tree.get_path(node)

        # Get all mutations relative to root
        muts = [branch.branch_attrs['mutations'].get(args.gene, []) for branch in path]
        # Flatten the list of mutations
        muts = [item for sublist in muts for item in sublist]
        # Get sequence at node
        node_seq = apply_muts_to_root(root_seq, muts)
        # Strip trailing stop codons
        stripped_seq = Seq(str(node_seq).rstrip('*'))
        # Strip hCoV-19/ from beginning of strain name
        strain = node.name.removeprefix('hCoV-19/')
        # Only keep records without stop codons (*)
        if '*' not in stripped_seq:
            sequence_records.append(SeqRecord(stripped_seq, strain, '', ''))

    # Write sequences to output FASTA file
    SeqIO.write(sequence_records, args.output, "fasta")
