#!/usr/bin/env python

import pandas as pd
import json
import argparse

def process_dataframe(df):
    # Separate metadata columns and standard columns
    standard_cols = [col for col in df.columns if not col.startswith('metadata_')]
    metadata_cols = [col for col in df.columns if col.startswith('metadata_')]

    # Create standard data dictionary
    data = df[standard_cols].to_dict(orient='records')

    # Process metadata if present
    if metadata_cols:
        metadata_data = df[metadata_cols].to_dict(orient='records')
        # Rename keys to remove 'metadata_' prefix and combine with main data
        for record, metadata in zip(data, metadata_data):
            cleaned_metadata = {k.replace('metadata_', ''): v for k, v in metadata.items()}
            record['metadata'] = cleaned_metadata

    return data

def read_excel_to_dict(filepath):
    data = {}
    xls = pd.ExcelFile(filepath)
    for sheet_name in xls.sheet_names:
        # Skip processing for specific sheets if needed
        if sheet_name.startswith('cv_'):
            continue

        df = xls.parse(sheet_name, skiprows=[1])

        # Remove columns that are not needed
        df = df.loc[:, ~df.columns.str.startswith('Unnamed')]

        df.fillna('', inplace=True)

        df.drop_duplicates(inplace=True)

        data[sheet_name] = process_dataframe(df)

    return data

def read_tsv_to_dict(filenames):
    data = {}
    for filename in filenames:
        entity_name = filename.split('.')[0]  # Assumes filename format "entity.tsv"
        df = pd.read_csv(filename, delimiter='\t', skiprows=[1])
        df.drop_duplicates(inplace=True)
        df.fillna('', inplace=True)
        data[entity_name] = process_dataframe(df)

    return data

def main(filepaths, output_json_path, is_excel=False):
    if is_excel:
        data = read_excel_to_dict(filepaths[0])  # Expects a list with one Excel file path
    else:
        data = read_tsv_to_dict(filepaths)  # Expects a list of TSV file paths

    with open(output_json_path, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert TSVs or an Excel file to JSON.")
    parser.add_argument("filepaths", nargs='+', help="Path(s) to the input TSV or single Excel file")
    parser.add_argument("output_json_path", help="Path to the output JSON file")
    parser.add_argument("--is_excel", action='store_true', help="Flag if the input is an Excel file")
    args = parser.parse_args()

    main(args.filepaths, args.output_json_path, args.is_excel)
