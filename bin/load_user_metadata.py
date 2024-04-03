#!/usr/bin/env python

import pandas as pd
import json
import argparse


def process_dataframe(df):
    # Separate metadata columns and standard columns
    standard_cols = [col for col in df.columns if not col.startswith("metadata_")]
    metadata_cols = [col for col in df.columns if col.startswith("metadata_")]

    # enforce str for taxon_id
    if "taxon_id" in standard_cols:
        df["taxon_id"] = df["taxon_id"].astype(str)

    # Create standard data dictionary
    data = df[standard_cols].to_dict(orient="records")

    # Process metadata if present
    if metadata_cols:
        metadata_data = df[metadata_cols].to_dict(orient="records")
        # Rename keys to remove 'metadata_' prefix and combine with main data
        for record, metadata in zip(data, metadata_data):
            cleaned_metadata = {
                k.replace("metadata_", ""): v for k, v in metadata.items()
            }
            record["metadata"] = cleaned_metadata

    return data


def read_excel_to_dict(filepath):
    data = {}
    xls = pd.ExcelFile(filepath)
    for sheet_name in xls.sheet_names:
        # Skip processing for specific sheets if needed
        if sheet_name.startswith("cv_"):
            continue

        df = xls.parse(sheet_name, skiprows=[1])

        # Remove columns that are not needed
        df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

        df.fillna("", inplace=True)

        df.drop_duplicates(inplace=True)

        data[sheet_name] = process_dataframe(df)

    return data


def read_tsv_to_dict(filenames):
    data = {}
    for filename in filenames:
        entity_name = filename.split(".")[0]  # Assumes filename format "entity.tsv"
        df = pd.read_csv(filename, delimiter="\t", skiprows=[1])
        df.drop_duplicates(inplace=True)
        df.fillna("", inplace=True)
        data[entity_name] = process_dataframe(df)

    return data


def extract_fastq_info(data, output_json_path):
    # Create a dictionary for fast access to experiments by their accession number
    experiments = {exp["accession"]: exp for exp in data["experiment"]}

    # Create a dictionary for fast access to samples by their accession number
    samples = {sample["accession"]: sample for sample in data["sample"]}

    output_data = []
    for run in data["run"]:
        experiment = experiments.get(run["experiment_accession"], {})
        sample = samples.get(experiment.get("sample_accession"), {})

        md5_1 = run["metadata"].get("fastq_md5", "").split(";")[0]
        if len(run["metadata"].get("fastq_md5", "").split(";")) > 1:
            md5_2 = run["metadata"].get("fastq_md5", "").split(";")[1]
        else:
            md5_2 = ""
        sample_dict = {
            "sample": run["experiment_accession"],
            "fastq_1": run.get("file_1", ""),
            "fastq_2": run.get("file_2", ""),
            "run_accession": run.get("accession", ""),
            "experiment_accession": run.get("experiment_accession", ""),
            "sample_accession": experiment.get("sample_accession", ""),
            "secondary_sample_accession": sample.get("secondary_accession", ""),
            "study_accession": experiment.get("study_accession", ""),
            "secondary_study_accession": experiment.get(
                "secondary_study_accession", ""
            ),
            "submission_accession": experiment.get("submission_accession", ""),
            "run_alias": run.get("alias", ""),
            "experiment_alias": experiment.get("alias", ""),
            "sample_alias": sample.get("alias", ""),
            "study_alias": experiment.get("study_alias", ""),
            "library_layout": experiment.get("library_layout", ""),
            "library_selection": experiment.get("library_selection", ""),
            "library_source": experiment.get("library_source", ""),
            "library_strategy": experiment.get("library_strategy", ""),
            "library_name": experiment.get("library_name", ""),
            "instrument_model": experiment.get("instrument_model", ""),
            "instrument_platform": experiment.get("platform", ""),
            "base_count": run.get("base_count", ""),
            "read_count": run.get("read_count", ""),
            "taxon_id": str(sample.get("tax_id", "")),  # Map tax_id to taxon_id
            "scientific_name": sample.get("scientific_name", ""),
            "sample_title": sample.get("title", ""),
            "experiment_title": experiment.get("title", ""),
            "study_title": experiment.get("study_title", ""),
            "sample_description": sample.get("description", ""),
            "fastq_md5": run.get("fastq_md5", ""),
            "fastq_bytes": run.get("fastq_bytes", ""),
            "fastq_ftp": run.get("fastq_ftp", ""),
            "fastq_galaxy": run.get("fastq_galaxy", ""),
            "fastq_aspera": run.get("fastq_aspera", ""),
            "cell_type": sample.get("cell_type", ""),
            "tissue_type": sample.get("tissue_type", ""),
            "cell_line": sample.get("cell_line", ""),
            "md5_1": md5_1,
            "md5_2": md5_2,
        }
        output_data.append(sample_dict)

    # Save the transformed data to a new JSON file
    with open(output_json_path, "w") as outfile:
        json.dump(output_data, outfile, indent=4)


def main(
    filepaths, output_metadata_json_path, output_samplesheet_json_path, is_excel=False
):
    if is_excel:
        data = read_excel_to_dict(
            filepaths[0]
        )  # Expects a list with one Excel file path
    else:
        data = read_tsv_to_dict(filepaths)  # Expects a list of TSV file paths

    with open(output_metadata_json_path, "w") as jsonfile:
        json.dump(data, jsonfile, indent=4)

    extract_fastq_info(data, output_samplesheet_json_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert TSVs or an Excel file to JSON."
    )
    parser.add_argument(
        "filepaths", nargs="+", help="Path(s) to the input TSV or single Excel file"
    )
    parser.add_argument(
        "output_metadata_json_path", help="Path to the output metadata JSON file"
    )
    parser.add_argument(
        "output_samplesheet_json_path", help="Path to the output samplesheet JSON file"
    )
    parser.add_argument(
        "--is_excel", action="store_true", help="Flag if the input is an Excel file"
    )
    args = parser.parse_args()

    main(
        args.filepaths,
        args.output_metadata_json_path,
        args.output_samplesheet_json_path,
        args.is_excel,
    )
