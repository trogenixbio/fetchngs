#!/usr/bin/env python

import json
import argparse

def update_run_info(input_json_path, output_json_path, cloud_prefix):
    # Load input JSON data
    with open(input_json_path, 'wr') as infile:
        data = json.load(infile)

    updated_run = []
    for run in data["run"]:
        # "file_1": "az://raw/results/fastq/DRX026011_DRR028935_1.fastq.gz",
        file_1 = run["file_1"].split("/")[-1]
        file_2 = run["file_1"].split("/")[-1]
        run["file_1"] = f"{cloud_prefix}/fastq/"
        run["file_2"] = f"{cloud_prefix}/fastq/"
        updated_run.append(run)

    # Update run information
    data["run"] = updated_run

    # Save the transformed data to a new JSON file
    with open(output_json_path, 'w') as outfile:
        json.dump(data, outfile, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract FASTQ file information from JSON based on experiment accession.")
    parser.add_argument("input_json_path", help="Path to the input JSON file")
    parser.add_argument("output_json_path", help="Path to the output JSON file")
    parser.add_argument("cloud_prefix_path", help="Path to the cloud prefix to update")
    args = parser.parse_args()

    update_run_info(args.input_json_path, args.output_json_path, args.cloud_prefix)
