#!/usr/bin/env python3

import argparse
import pandas as pd

def annotate_sequence_counts(input_file, output_file, dataset_suffix):
    # Read the TSV input file
    df = pd.read_csv(input_file, sep='\t')

    # Append dataset_suffix to the 'location' column
    df['location'] = df['location'] + f"_{dataset_suffix}"

    # Modify 'variant' where the value is 'other'
    df.loc[df['variant'] == 'other', 'variant'] = df['variant'] + f"_{dataset_suffix}"

    # Write the annotated DataFrame to the output file
    df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Annotate sequence counts with dataset suffix.")
    parser.add_argument('--input', required=True, help='Path to the input TSV file')
    parser.add_argument('--output', required=True, help='Path to the output TSV file')
    parser.add_argument('--dataset-suffix', required=True, help='Suffix to append to location and variant (if "other")')

    args = parser.parse_args()

    annotate_sequence_counts(args.input, args.output, args.dataset_suffix)

if __name__ == "__main__":
    main()
