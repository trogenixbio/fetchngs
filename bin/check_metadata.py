#!/usr/bin/env python

import json
import re
import argparse
import sys


def load_json(filepath):
    with open(filepath, "r") as file:
        return json.load(file)


def validate_output_against_schema(output_data, schema):
    validation_report = []
    validation_errors = []

    for entity, records in output_data.items():
        if entity not in schema:
            error_msg = f"Entity '{entity}' not found in schema."
            validation_errors.append(error_msg)
            validation_report.append(error_msg)
            continue

        entity_schema = schema[entity]["fields"]
        schema_fields = {field["name"]: field for field in entity_schema}

        for record in records:
            for field_name, field_info in schema_fields.items():
                value = record.get(field_name)
                if field_info["cardinality"] == "mandatory" and value is None:
                    error_msg = f"FAILURE mandatory field '{field_name}' missing in entity '{entity}'."
                    validation_errors.append(error_msg)
                elif value is not None:
                    if field_info["field_type"] == "TEXT_FIELD" and not isinstance(
                        value, str
                    ):
                        error_msg = f"FAILURE field '{field_name}' in entity '{entity}' should be a string."
                        validation_errors.append(error_msg)

                    if (
                        field_info["field_type"] == "TEXT_CHOICE_FIELD"
                        and value not in field_info["cv"]
                    ):
                        if not (
                            field_info["cardinality"] == "optional"
                            and (value is None or value == "")
                        ):
                            error_msg = f"FAILURE field '{field_name}' in entity '{entity}' has value '{value}' not in controlled vocabulary {field_info['cv']}."
                            validation_errors.append(error_msg)

                    if field_info.get("regex"):
                        if not (
                            field_info["cardinality"] == "optional"
                            and (value is None or value == "")
                        ):
                            if not re.match(field_info["regex"], value):
                                error_msg = f"Field '{field_name}' in entity '{entity}' with value '{value}' does not match regex pattern {field_info['regex']}."
                                validation_errors.append(error_msg)

                # Success report for each field
                validation_report.append(
                    f"SUCCESS field '{field_name}' in entity '{entity}' passes validation."
                )

    return "\n".join(validation_report), validation_errors


def save_valid_metadata(output_data):
    with open("metadata-valid.json", "w") as file:
        json.dump(output_data, file, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate JSON output against a schema."
    )
    parser.add_argument("output_json", help="Path to the output JSON file")
    parser.add_argument("schema_json", help="Path to the JSON schema file")

    args = parser.parse_args()

    schema = load_json(args.schema_json)
    output_data = load_json(args.output_json)

    report, errors = validate_output_against_schema(output_data, schema)

    if errors:
        print("\nValidation Errors Found:")
        for error in errors:
            print(error)
        sys.exit(
            f"Validation failed with {len(errors)} errors. Please see validation_report.txt for details."
        )
    else:
        print("Validation successful. No errors found.")
        save_valid_metadata(output_data)
        print("The validated metadata has been saved to metadata-valid.json.")
        print(report)
