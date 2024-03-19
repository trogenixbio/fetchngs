#!/usr/bin/env python

import csv
import json
import sys

def tsv_to_json(tsv_filepath, json_filepath):
    with open(tsv_filepath, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter=',')
        data = list(reader)

    with open(json_filepath, 'w') as jsonfile:
        json.dump(data, jsonfile, ensure_ascii=False, indent=4)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python sample_to_json.py <input_path> <output_path>")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_json = sys.argv[2]

    tsv_to_json(input_tsv, output_json)
    print(f"Converted {input_tsv} to {output_json} successfully.")
