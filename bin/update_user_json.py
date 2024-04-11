#!/usr/bin/env python

import json
import argparse


def update_run_info_metadata(
    input_metadata_json_path, output_metadata_json_path, cloud_prefix, pub_internal
):
    # Load input JSON data
    with open(input_metadata_json_path, "r") as infile:
        data = json.load(infile)

    # Create a dictionary for fast access to experiments by their accession number
    experiments = {exp["accession"]: exp for exp in data["experiment"]}

    updated_run = []
    for run in data["run"]:
        experiment = experiments.get(run["experiment_accession"], {})
        run["file_1"] = (
            f"{cloud_prefix}/{pub_internal}/{experiment['study_accession']}/{experiment['library_strategy']}/fastq/{run['experiment_accession']}_{run['accession']}.fastq.gz"
        )
        run["file_2"] = (
            f"{cloud_prefix}/{pub_internal}/{experiment['study_accession']}/{experiment['library_strategy']}/fastq/{run['experiment_accession']}_{run['accession']}.fastq.gz"
        )
        updated_run.append(run)

    # Update run information
    data["run"] = updated_run

    # Save the transformed data to a new JSON file
    with open(output_metadata_json_path, "w") as outfile:
        json.dump(data, outfile, indent=4)


def update_run_info_samplesheet(
    input_samplesheet_json_path,
    output_samplesheet_json_path,
    cloud_prefix,
    pub_internal,
):
    # Load the JSON data
    with open(input_samplesheet_json_path, "r") as file:
        data = json.load(file)

    # Modify the fastq_1 and fastq_2 paths
    for item in data:
        if "fastq_1" in item:
            fastq_1 = item["fastq_1"].split("/")[-1]
            item["fastq_1"] = (
                f"{cloud_prefix}/{pub_internal}/{item['study_accession']}/{item['library_strategy']}/fastq/{item['experiment_accession']}_{item['run_accession']}.fastq.gz"
            )
        if "fastq_2" in item:
            fastq_2 = item["fastq_2"].split("/")[-1]
            item["fastq_2"] = (
                f"{cloud_prefix}/{pub_internal}/{item['study_accession']}/{item['library_strategy']}/fastq/{item['experiment_accession']}_{item['run_accession']}.fastq.gz"
            )

    # Save the modified data to a new JSON file
    with open(output_samplesheet_json_path, "w") as file:
        json.dump(data, file, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract FASTQ file information from JSON based on experiment accession."
    )
    parser.add_argument(
        "input_metadata_json_path", help="Path to the input metadata JSON file"
    )
    parser.add_argument(
        "input_samplesheet_json_path", help="Path to the input samplesheet JSON file"
    )
    parser.add_argument(
        "output_metadata_json_path", help="Path to the input metadata JSON file"
    )
    parser.add_argument(
        "output_samplesheet_json_path", help="Path to the output samplesheet JSON file"
    )
    parser.add_argument("cloud_prefix", help="Path to the cloud prefix to update")
    parser.add_argument("pub_internal", help="Internal/public out dir")

    args = parser.parse_args()

    update_run_info_metadata(
        args.input_metadata_json_path,
        args.output_metadata_json_path,
        args.cloud_prefix,
        args.pub_internal,
    )
    update_run_info_samplesheet(
        args.input_samplesheet_json_path,
        args.output_samplesheet_json_path,
        args.cloud_prefix,
        args.pub_internal,
    )
