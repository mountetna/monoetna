"""
Ideally the output JSON data files are in the same format as the Timur
consignments, which means we can re-use UI components.
This format converts a data frame / vector from


        col1    col2    col3
row1     1        3      5.1
row2     0       -8      87
row3     99      1       12.3

to an object like:

{
    "col_names": ["col1", "col2", "col3"],
    "num_cols": 3,
    "num_rows": 3,
    "row_names": ["row1", "row2", "row3"],
    "rows": [
        [1, 3, 5.1],
        [0, -8, 87],
        [99, 1, 12.3]
    ]
}
"""

import argparse
import csv
import json
import tempfile


def convert_tsv_to_json(input, output):
    with open(input, "r") as input_file, open(output, "w") as output_file:
        reader = csv.reader(input_file, delimiter="\t")
        output_data = {
            "rows": []
        }
        headers = False
        row_names = []
        for row in reader:
            row_data = row[1:]
            if not headers:
                headers = True
                output_data["col_names"] = row_data
                output_data["num_cols"] = len(row_data)
                continue
            row_names.append(row[0])

            numeric_row = []
            for datum in row_data:
                try:
                    numeric_row.append(float(datum))
                except:
                    numeric_row.append(datum)

            output_data["rows"].append(numeric_row)
        output_data["num_rows"] = len(row_names)
        output_data["row_names"] = row_names
        json.dump(output_data, output_file, indent=2)
        return output_data


def test_convert():
    results = convert_tsv_to_json(
        "fixtures/data.tsv", tempfile.NamedTemporaryFile().name)
    assert (results == {
        "num_rows": 2,
        "row_names": [
            "row1",
            "row2"
        ],
        "num_cols": 2,
        "rows": [
            [
                1, 2
            ],
            [
                2, 1
            ]
        ],
        "col_names": [
            "col1",
            "col2"
        ]
    })


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert TSV to Consignment JSON format')
    parser.add_argument('--input', '-i', type=str,
                        default='data.tsv',
                        help='Input TSV file',
                        dest='input')
    parser.add_argument('--output', '-o', type=str,
                        default='data.json',
                        help='Output JSON file',
                        dest='output')

    arguments = parser.parse_args()
    convert_tsv_to_json(arguments.input, arguments.output)
