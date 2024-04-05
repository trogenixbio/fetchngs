#!/usr/bin/env python

import json
import argparse
import csv


def load_mappings(mapping_json_path):
    with open(mapping_json_path, "r") as file:
        mappings = json.load(file)
    return mappings["field_mapping"], mappings["group_mapping"]


def apply_field_mapping(item, schema, field_mapping, group_mapping):
    mapped_item = {entity: {} for entity in schema.keys()}
    fieldnames_schema = {entity: {} for entity in schema.keys()}

    # First, apply direct mappings based on the schema field names
    for entity in schema:
        for field_info in schema[entity]["fields"]:
            field_name = field_info["name"]
            if "fields" in fieldnames_schema[entity]:
                fieldnames_schema[entity]["fields"].append(field_name)
            else:
                fieldnames_schema[entity]["fields"] = [field_name]
            if field_name in item:
                mapped_item[entity][field_name] = item[field_name]

    # Then, apply explicit field mappings
    for input_field, value in item.items():
        if input_field in field_mapping:
            if not any(
                isinstance(el, list) for el in field_mapping[input_field]
            ):  # If a singlet
                entity, schema_field = field_mapping[input_field]
                if schema_field not in fieldnames_schema[entity]["fields"]:
                    if "metadata" not in mapped_item[entity]:
                        mapped_item[entity]["metadata"] = {}
                    mapped_item[entity]["metadata"][schema_field] = value
                else:
                    mapped_item[entity][schema_field] = value
            else:  # If a multimapper
                for mapping in field_mapping[input_field]:
                    entity, schema_field = mapping
                    # Create metadata sub-dictionary if a direct match is not found
                    if schema_field not in fieldnames_schema[entity]["fields"]:
                        if "metadata" not in mapped_item[entity]:
                            mapped_item[entity]["metadata"] = {}
                        mapped_item[entity]["metadata"][schema_field] = value
                    else:
                        mapped_item[entity][schema_field] = value
        elif input_field in group_mapping:
            entity = group_mapping[input_field]
            if "metadata" not in mapped_item[entity]:
                mapped_item[entity]["metadata"] = {}
            mapped_item[entity]["metadata"][input_field] = value
        else:
            print(f"Field missing from mapping {input_field}")

    return mapped_item


def convert_to_schema_format(input_data, schema, field_mapping, group_mapping):
    output = {key: [] for key in schema.keys()}
    for item in input_data:
        mapped_item = apply_field_mapping(item, schema, field_mapping, group_mapping)

        for entity, entity_def in schema.items():
            entity_obj = {}
            # Ensure all schema-defined fields are accounted for, even if not present in the input
            for field_info in entity_def["fields"]:
                field_name = field_info["name"]
                # If the field is present in the mapped item, use that value
                if field_name in mapped_item[entity]:
                    entity_obj[field_name] = mapped_item[entity][field_name]
                # If the field is missing, use an empty string as a placeholder
                else:
                    entity_obj[field_name] = ""

            # Include any additional metadata fields that were added during mapping
            if "metadata" in mapped_item[entity]:
                # Ensure metadata is not overwritten by empty strings if it's already populated
                if "metadata" not in entity_obj or not entity_obj["metadata"]:
                    entity_obj["metadata"] = mapped_item[entity]["metadata"]

            # Append the entity object if it has any fields defined
            if entity_obj:
                output[entity].append(entity_obj)

    return output


def json_to_tsv(json_data, tsv_filepath):
    if not json_data:
        print("No data to write to TSV.")
        return

    written_items = []  # List to track items already processed
    with open(tsv_filepath, "w", newline="") as tsvfile:
        writer = None
        for item in json_data:
            flat_item = {k: v for k, v in item.items() if k != "metadata"}
            if "metadata" in item:
                for meta_key, meta_value in item["metadata"].items():
                    flat_item[f"metadata_{meta_key}"] = meta_value

            # Convert the dictionary to a tuple of items for comparison
            item_tuple = tuple(flat_item.items())
            if item_tuple not in written_items:
                written_items.append(item_tuple)  # Track this item as written
                if writer is None:
                    headers = flat_item.keys()
                    writer = csv.DictWriter(tsvfile, fieldnames=headers, delimiter="\t")
                    writer.writeheader()
                writer.writerow(flat_item)

    print(f"Data successfully written to {tsv_filepath}")


def main(input_json_path, schema_json_path, output_json_path, mapping_json_path):

    input_data = json.load(open(input_json_path, "r"))
    schema = json.load(open(schema_json_path, "r"))

    field_mapping, group_mapping = load_mappings(mapping_json_path)

    output_data = convert_to_schema_format(
        input_data, schema, field_mapping, group_mapping
    )
    with open(output_json_path, "w") as file:
        json.dump(output_data, file, indent=4, ensure_ascii=False)

    for entity, records in output_data.items():
        # Define a TSV file path for each entity
        tsv_filepath = f"{entity}.tsv"

        # Flatten the records if necessary
        # This step depends on the structure of your records
        # For now, assuming records is a list of flat dictionaries
        # flat_records = records  # Implement flattening if necessary

        # # Convert the JSON data to TSV
        # json_to_tsv(flat_records, tsv_filepath)
        json_to_tsv(records, tsv_filepath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert JSON to conform with a given schema using external field and group mapping."
    )
    parser.add_argument("input_json", help="Path to the input JSON file")
    parser.add_argument("schema_json", help="Path to the JSON schema file")
    parser.add_argument("output_json", help="Path to the output JSON file")
    parser.add_argument("mapping_json", help="Path to the JSON mapping file")

    args = parser.parse_args()

    main(args.input_json, args.schema_json, args.output_json, args.mapping_json)
