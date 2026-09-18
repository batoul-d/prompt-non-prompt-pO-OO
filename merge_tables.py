#!/usr/bin/env python3

import csv
import argparse


def is_empty(row):
    """Return True if a row is completely empty."""
    return not row or all(cell.strip() == "" for cell in row)


parser = argparse.ArgumentParser(
    description="Merge multiple CSV table files, keeping one header and all configuration names."
)

parser.add_argument(
    "-o",
    "--output",
    default="merged_table.csv",
    help="Output file name (default: merged_table.csv)"
)

parser.add_argument(
    "files",
    nargs="+",
    help="Input table files to merge"
)

args = parser.parse_args()

header_written = False

with open(args.output, "w", newline="") as outfile:
    writer = csv.writer(
        outfile,
        delimiter=",",
        quotechar='"',
        quoting=csv.QUOTE_MINIMAL,
        lineterminator="\n"
    )

    for filename in args.files:

        print(f"Adding: {filename}")

        with open(filename, "r", newline="") as infile:
            reader = csv.reader(
                infile,
                delimiter=",",
                quotechar='"'
            )

            first_nonempty_row = True

            for row in reader:

                # Ignore completely empty rows
                if is_empty(row):
                    continue

                # First non-empty row = header
                if not header_written:
                    writer.writerow(row)
                    header_written = True
                    first_nonempty_row = False
                    continue

                # For every subsequent file, skip its first non-empty row
                # (the header)
                if first_nonempty_row:
                    first_nonempty_row = False
                    continue

                # Write everything else
                writer.writerow(row)

print(f"\nMerged table written to: {args.output}")