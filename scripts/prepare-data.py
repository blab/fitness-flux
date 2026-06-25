"""
Prepare case counts and clade counts data for analysis by subsetting to
specific date range and locations, pruning recent clade counts, and collapsing clades.
"""
import argparse
import pandas as pd
import sys

from datetime import datetime, timedelta

SEQ_COUNTS_DTYPES = {
    'location': 'string',
    'clade': 'string',
    'sequences': 'int64',
}

# Default min date is 1 year before today's date
DEFAULT_MIN_DATE = (datetime.now() - timedelta(days=365)).strftime('%Y-%m-%d')

# Default max date is today's date
DEFAULT_MAX_DATE = (datetime.now() - timedelta(days=1)).strftime('%Y-%m-%d')

def positive_int(value):
    """
    Custom argparse type function to verify only
    positive integers are provided as arguments
    """
    int_value = int(value)
    if int_value <= 0:
        print(f"ERROR: {int_value} is not a positive integer.", file=sys.stderr)
        sys.exit(1)
    return int_value


def unit_frequency(value):
    """
    Custom argparse type function to verify a frequency strictly between 0 and 1.
    """
    float_value = float(value)
    if not (0.0 < float_value < 1.0):
        print(f"ERROR: {float_value} is not a frequency in the open interval (0, 1).", file=sys.stderr)
        sys.exit(1)
    return float_value


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--seq-counts", metavar="TSV", required=True,
        help="Path to clade counts TSV with four columns: 'location','clade','date','sequences'")
    parser.add_argument("--min-date", default=DEFAULT_MIN_DATE,
        help="The minimum cutoff for date (inclusive), formatted as 'YYYY-MM-DD'.\n"
             "(default: %(default)s)")
    parser.add_argument("--max-date", default=DEFAULT_MAX_DATE,
        help="The maximum cutoff for date (inclusive), formatted as 'YYYY-MM-DD'.\n"
             "(default: %(default)s)")
    parser.add_argument("--location-min-seq", type=positive_int, default=1,
        help="The mininum number of sequences a location must have within the "
             "days-min-seq to be included in analysis.\n"
             "(default: %(default)s)")
    parser.add_argument("--location-min-seq-days", type=positive_int,
        help="The number of days (counting back from the cutoff date) to use as the date range "
             "for counting the number of sequences per location to determine if a location is included in analysis.\n"
             "If not provided, will count sequences from all dates included in analysis date range.")
    parser.add_argument("--excluded-locations",
        help="File with a list locations to exclude from analysis.")
    parser.add_argument("--clade-min-count", type=positive_int,
        help="The minimum number of sequences a clade must have to be included as its own variant.\n"
             "All clades with fewer than the minimum will be collapsed as 'other'.")
    parser.add_argument("--clade-min-freq", type=unit_frequency,
        help="The minimum mean frequency (clade sequence count / window total, in the open\n"
             "interval (0, 1)) a clade must have to be included as its own variant.\n"
             "All clades below the minimum will be collapsed as 'other'. This is the\n"
             "frequency-based analogue of --clade-min-count that auto-scales with window size;\n"
             "if both are given a clade must satisfy both thresholds.")
    parser.add_argument("--clade-min-count-days", type=positive_int,
        help="The number of days (counting back from the cutoff date) to use as the date range "
             "for counting the number of sequences per clade to determine if a clade is included as its own variant.\n"
             "If not provided, will count sequences from all dates included in analysis date range.")
    parser.add_argument("--force-include-clades", nargs="*",
        help="Clades to force include in the output regardless of sequences counts.")
    parser.add_argument("--force-exclude-clades", nargs="*",
        help="Clades to force exclude in the output regardless of sequences counts.")
    parser.add_argument("--output-seq-counts", required=True,
        help="Path to output TSV file for the prepared variants data.")

    args = parser.parse_args()

    ###########################################################################
    ##################### Rules for subsetting by date ########################
    ###########################################################################
    # Max date to include for analysis
    print(f"Setting max date (inclusive) as {args.max_date!r}.")
    max_date = datetime.strptime(args.max_date, '%Y-%m-%d')

    print(f"Setting min date (inclusive) as {args.min_date!r}.")
    min_date = datetime.strptime(args.min_date, '%Y-%m-%d')

    ###########################################################################
    ################### Rules for subsetting by location ######################
    ###########################################################################
    # Load entire clade counts data since we need to to find all locations
    seq_counts = pd.read_csv(args.seq_counts, sep='\t', parse_dates=['date'], dtype=SEQ_COUNTS_DTYPES)

    # Set the min_date as the default min date for counting sequences per location
    # to count sequences per location over the entire analysis date range
    min_location_seq_date = min_date
    if args.location_min_seq_days is not None:
        print(
            f"Only including locations that have at least {args.location_min_seq} sequence(s)",
            f"in the last {args.location_min_seq_days} days of the analysis date range."
        )
        # Calculate the minimum date for sequences per location as *location_min_seq_days* before the max date
        # Subtract 1 from days in calculation since we are including the max date
        min_location_seq_date = max_date - timedelta(days=(args.location_min_seq_days - 1))
    else:
        print(
            f"Only including locations that have at least {args.location_min_seq} sequence(s)",
            "in the analysis date range."
        )

    # Subset to locations and sequences within the date range from min_location_seq_date to max_date
    # and group by location to get the total number of sequences per location in this date range
    location_and_seqs = seq_counts.loc[
        (seq_counts['date'] >= min_location_seq_date if min_location_seq_date else True) &
        (seq_counts['date'] <= max_date),
        ['location', 'sequences']
    ]
    seqs_per_location = location_and_seqs.groupby(['location'], as_index=False).sum()

    # Get a set of locations that meet the location_min_seq requirement
    locations_with_min_seq = set(seqs_per_location.loc[seqs_per_location['sequences'] >= args.location_min_seq, 'location'])

    # Load manually annotated excluded locations if provided
    excluded_locations = set()
    if args.excluded_locations:
        with open(args.excluded_locations, 'r') as f:
            excluded_locations = {line.rstrip() for line in f}

        print(f"Excluding the following requested locations: {sorted(excluded_locations)}.")

    # Remove excluded-locations from the set of locations to include in analysis
    locations_to_include = locations_with_min_seq - excluded_locations
    print(f"Locations that will be included: {sorted(locations_to_include)}.")

    assert len(locations_to_include) > 0, \
        "All locations have been excluded. Try again with different options, e.g. lowering the `--location-min-seq` cutoff.\n" + \
        f"Here's a summary of available sequences per location:\n{seqs_per_location.to_dict(orient='records')}"

    ###########################################################################
    ############## Rules for collapsing clades to variants ####################
    ###########################################################################

    # Duplicate 'clade' column to 'variant' column to faciliate mapping to 'other'
    seq_counts['variant'] = seq_counts.loc[:, 'clade']

    # Keep track of clades that are force included so that they can bypass the sequence counts check
    force_included_clades = set()
    if args.force_include_clades:
        for variant in args.force_include_clades:
            force_included_clades.add(variant)

        print(f"Force including the following variants: {args.force_include_clades}")

    # Keep track of clades that are force excluded so that they can bypass the sequence counts check
    force_excluded_clades = set()
    if args.force_exclude_clades:
        for variant in args.force_exclude_clades:
            force_excluded_clades.add(variant)

        print(f"Force excluding the following variants: {args.force_exclude_clades}")

    # Collapse small clades into "other" if a count and/or frequency threshold is provided
    if args.clade_min_count or args.clade_min_freq:
        # Set the min_date as the default min date for counting sequences per clade
        # to count sequences per clade over the entire analysis date range
        min_clades_seq_date = min_date
        if args.clade_min_count_days is not None:
            window_description = f"in the last {args.clade_min_count_days} days of the analysis date range"
            # Calculate the minimum date for sequences per clade as *clade_min_count_days* before the max date
            # Subtract 1 from days in calculation since we are including the max date
            min_clades_seq_date = max_date - timedelta(days=(args.clade_min_count_days - 1))
        else:
            window_description = "in the analysis date range (inclusive)"

        # Subset to clades and sequences within the date range from min_clades_seq_date to max_date
        # and group by clade to get the total number of sequences per clade in this date range
        clades_and_seqs = seq_counts.loc[
            (seq_counts['date'] >= min_clades_seq_date if min_clades_seq_date else True) &
            (seq_counts['date'] <= max_date),
            ['clade', 'sequences']
        ]
        seqs_per_clade = clades_and_seqs.groupby(['clade'], as_index=False).sum()
        total_seqs = seqs_per_clade['sequences'].sum()

        # A clade is kept only if it satisfies every provided threshold
        clades_to_keep = set(seqs_per_clade['clade'])
        if args.clade_min_count:
            print(
                f"Collapsing clades that have less than {args.clade_min_count} sequence(s)",
                f"{window_description} into a single 'other' variant."
            )
            clades_to_keep &= set(seqs_per_clade.loc[seqs_per_clade['sequences'] >= args.clade_min_count, 'clade'])
        if args.clade_min_freq:
            # Mean frequency over the window; the equivalent count auto-scales with window size
            min_freq_seqs = args.clade_min_freq * total_seqs
            print(
                f"Collapsing clades with a mean frequency below {args.clade_min_freq}",
                f"(fewer than {min_freq_seqs:.1f} of {total_seqs} sequences) {window_description} into a single 'other' variant."
            )
            clades_to_keep &= set(seqs_per_clade.loc[seqs_per_clade['sequences'] >= min_freq_seqs, 'clade'])

        # Remove force excluded clades from the keep list
        clades_to_keep -= force_excluded_clades

        # Replace variant with 'other' if they are not force included and do not meet the threshold(s)
        seq_counts.loc[~seq_counts['clade'].isin(force_included_clades | clades_to_keep), 'variant'] = 'other'

    # Replace 'recombinant' clade with 'other'
    seq_counts.loc[seq_counts['clade'].isin(['recombinant']), 'variant'] = 'other'

    # Add clades with unknown variants to 'other' variant group
    seq_counts.loc[pd.isna(seq_counts['variant']), 'variant'] = 'other'

    # Collapse the variants of the same location and date
    seq_counts = seq_counts.groupby(['location', 'variant', 'date'], as_index=False).sum(numeric_only=False)
    seq_counts.drop(columns=['clade'], inplace=True)


    ###########################################################################
    ##################### Rules for pruning sequence counts ###################
    ###########################################################################

    # The default max date for clade counts is the max date
    max_clade_date = max_date

    ###########################################################################
    ########################## Subset and output data #########################
    ###########################################################################

    # Subset the clade counts data by date and locations
    seq_counts = seq_counts.loc[
        (seq_counts['date'] >= min_date if min_date else True) &
        (seq_counts['date'] <= max_clade_date) &
        (seq_counts['location'].isin(locations_to_include))
    ]

    included_variants = seq_counts['variant'].unique()
    print(f"Variants that will be included: {sorted(included_variants)}.")

    assert len(included_variants) > 0, \
        "All variants have been excluded. Try again with different options, e.g. lowering the `--clade-min-count` cutoff."

    # Sort variants subset and print to output file
    seq_counts.sort_values(['location', 'variant', 'date']) \
          .to_csv(
              args.output_seq_counts,
              sep='\t',
              index=False,
          )
