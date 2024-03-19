#!/usr/bin/env python

import json
import argparse

def load_mappings(mapping_json_path):
    with open(mapping_json_path, 'r') as file:
        mappings = json.load(file)
    return mappings['field_mapping'], mappings['group_mapping']

def apply_field_mapping(item, schema, field_mapping, group_mapping):
    mapped_item = {entity: {} for entity in schema.keys()}
    for input_field, value in item.items():
        if input_field in field_mapping:
            entity, schema_field = field_mapping[input_field]
            mapped_item[entity][schema_field] = value
        elif input_field in group_mapping:
            entity = group_mapping[input_field]
            if "metadata" not in mapped_item[entity]:
                mapped_item[entity]["metadata"] = {}
            mapped_item[entity]["metadata"][input_field] = value
    return mapped_item

def convert_to_schema_format(input_data, schema, field_mapping, group_mapping):
    output = {key: [] for key in schema.keys()}
    for item in input_data:
        mapped_item = apply_field_mapping(item, schema, field_mapping, group_mapping)
        for entity, entity_fields in mapped_item.items():
            entity_obj = {}
            for field in schema[entity]["fields"]:
                field_name = field["name"]
                if field_name in entity_fields:
                    entity_obj[field_name] = entity_fields[field_name]
                elif field["cardinality"] == "mandatory":
                    entity_obj[field_name] = "not provided"
            if "metadata" in entity_fields:
                entity_obj["metadata"] = entity_fields["metadata"]
            output[entity].append(entity_obj)
    return output

def main(input_json_path, schema_json_path, output_json_path, mapping_json_path):
    input_data = json.load(open(input_json_path, 'r'))
    schema = json.load(open(schema_json_path, 'r'))
    field_mapping, group_mapping = load_mappings(mapping_json_path)

    output_data = convert_to_schema_format(input_data, schema, field_mapping, group_mapping)
    with open(output_json_path, 'w') as file:
        json.dump(output_data, file, indent=4, ensure_ascii=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert JSON to conform with a given schema using external field and group mapping.")
    parser.add_argument("input_json", help="Path to the input JSON file")
    parser.add_argument("schema_json", help="Path to the JSON schema file")
    parser.add_argument("output_json", help="Path to the output JSON file")
    parser.add_argument("mapping_json", help="Path to the JSON mapping file")

    args = parser.parse_args()

    main(args.input_json, args.schema_json, args.output_json, args.mapping_json)

# Define a mapping from input JSON fields to schema groups
# group_mapping = {
#     "sample_related_field": "sample",
#     "study_related_field": "study",
#     "experiment_related_field": "experiment",
#     "run_related_field": "run",
#     # Add other mappings as needed. These fields don't map to specific schema fields but to a group
# }

# # Define specific field mappings as before
# field_mapping = {
#     "sample_accession": ("sample", "alias"),
#     "study_accession": ("study", "accession"),
#     "study_title": ("study", "title"),
#     "experiment_accession": ("experiment", "accession"),
#     "experiment_title": ("experiment", "title"),
#     "run_accession": ("run", "accession"),
#     # Continue with specific field mappings
# }
