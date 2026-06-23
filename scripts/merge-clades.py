#!/usr/bin/env python3
"""
Relabel (merge) clade/variant values in a TSV read from stdin and written to
stdout. Used by the provisioning step to collapse several source clades into one
target clade (e.g. SARS-CoV-2 19A + 19B -> WT) before the rest of the workflow
keys off the relabeled column.

The merges are passed as a JSON object via --merges, of the form

    {column: {target: [source, source, ...], ...}, ...}

For each named column, any cell whose value is one of the sources is rewritten to
the target. An empty object is an identity pass-through (so a virus without a
merge configured is unaffected).
"""
import argparse
import csv
import json
import sys


def build_remap(merges):
    # column -> {source_value: target_value}
    return {
        column: {src: target for target, sources in groups.items() for src in sources}
        for column, groups in merges.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merges", required=True, help="JSON {column: {target: [sources]}}")
    args = parser.parse_args()

    remap = build_remap(json.loads(args.merges))

    reader = csv.DictReader(sys.stdin, delimiter="\t")
    writer = csv.DictWriter(
        sys.stdout, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in reader:
        for column, mapping in remap.items():
            if column in row and row[column] in mapping:
                row[column] = mapping[row[column]]
        writer.writerow(row)


if __name__ == "__main__":
    main()
