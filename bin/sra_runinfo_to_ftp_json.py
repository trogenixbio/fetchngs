#!/usr/bin/env python


import argparse
import json
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


def parse_args(args=None):
    Description = "Process JSON data to create a samplesheet with FTP download links and md5sums in JSON format."
    Epilog = "Example usage: python script.py input.json output.json"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in", help="Input JSON file with sample data.")
    parser.add_argument(
        "file_out",
        type=Path,
        help="Output JSON file with selected data and additional fields.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="Log level (default: WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(args)


def parse_json_input(file_in):
    with open(file_in, "r") as fin:
        data = json.load(fin)

    if not data:
        return []

    runinfo = []

    for row in data:
        experiment_accession = row.get("experiment_accession")
        if (
            not experiment_accession
        ):  # Check if 'experiment_accession' is missing or empty
            raise ValueError("Missing 'experiment_accession' in one or more items.")

        run_accession = row.get("run_accession")
        row_id = (
            f"{experiment_accession}_{run_accession}"
            if run_accession
            else experiment_accession
        )
        row["id"] = row_id

        fq_files = row["fastq_ftp"].split(";") if row.get("fastq_ftp") else []
        fq_md5 = row["fastq_md5"].split(";") if row.get("fastq_md5") else []
        if row.get("library_layout") == "SINGLE":
            row.update(
                {
                    "fastq_1": fq_files[0] if fq_files else None,
                    "fastq_2": None,
                    "md5_1": fq_md5[0] if fq_md5 else None,
                    "md5_2": None,
                    "single_end": "true",
                }
            )
        elif row.get("library_layout") == "PAIRED":
            row.update(
                {
                    "fastq_1": fq_files[0] if len(fq_files) > 0 else None,
                    "fastq_2": fq_files[1] if len(fq_files) > 1 else None,
                    "md5_1": fq_md5[0] if len(fq_md5) > 0 else None,
                    "md5_2": fq_md5[1] if len(fq_md5) > 1 else None,
                    "single_end": "false",
                }
            )

        runinfo.append(row)

    return runinfo


def write_output(runinfo, file_out):
    with open(file_out, "w", encoding="utf-8") as fout:
        json.dump(runinfo, fout, ensure_ascii=False, indent=4)


def main(args=None):
    args = parse_args(args)
    logging.basicConfig(
        level=args.log_level.upper(), format="[%(levelname)s] %(message)s"
    )

    file_in = Path(args.file_in)
    file_out = args.file_out

    if not file_in.is_file():
        logger.critical(f"Input file {file_in} not found.")
        sys.exit(1)

    file_out.parent.mkdir(parents=True, exist_ok=True)

    runinfo = parse_json_input(file_in)
    write_output(runinfo, file_out)


if __name__ == "__main__":
    sys.exit(main())
