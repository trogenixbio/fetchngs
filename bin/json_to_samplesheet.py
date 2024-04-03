#!/usr/bin/env python

import json
import csv
import argparse


def json_to_csv(input_path, output_path, include_cols=[], exclude_cols=[]):
    with open(input_path, "r") as json_file:
        data = json.load(json_file)

    with open(output_path, "w", newline="") as csv_file:
        if data:
            # Base columns to include if present in the JSON data
            basic_cols = ["sample", "fastq_1", "fastq_2"]
            # Determine all possible headers from the data
            all_headers = list(data[0].keys())

            # Start with basic columns, add include_cols next, and then append other columns not excluded
            ordered_headers = basic_cols + [
                col
                for col in include_cols
                if col in all_headers and col not in basic_cols
            ]
            ordered_headers += [
                col
                for col in all_headers
                if col not in ordered_headers and col not in exclude_cols
            ]

            writer = csv.DictWriter(
                csv_file, fieldnames=ordered_headers, quoting=csv.QUOTE_ALL
            )
            writer.writeheader()
            for item in data:
                # Only include the keys that are in our ordered list of headers
                filtered_item = {k: v for k, v in item.items() if k in ordered_headers}
                writer.writerow(filtered_item)


def main():
    parser = argparse.ArgumentParser(
        description="Convert JSON to CSV with inclusion and exclusion of specific columns."
    )
    parser.add_argument("input_json", help="Path to the input JSON file")
    parser.add_argument("output_csv", help="Path to the output CSV file")
    parser.add_argument(
        "--include",
        dest="pipeline_type",
        nargs="?",
        const="",
        default=None,
        help="Comma-separated column names to forcefully include after basic columns",
    )
    parser.add_argument(
        "--exclude",
        dest="exclude_cols",
        nargs="?",
        const="",
        default="",
        help="Comma-separated column names to exclude",
    )

    args = parser.parse_args()

    include_cols = args.pipeline_type.split(",") if args.pipeline_type else []
    exclude_cols = args.exclude_cols.split(",") if args.exclude_cols else []

    json_to_csv(args.input_json, args.output_csv, include_cols, exclude_cols)
    print(f"Modified JSON data has been saved to {args.output_csv}")


if __name__ == "__main__":
    main()
