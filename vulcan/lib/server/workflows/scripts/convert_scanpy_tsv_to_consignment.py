"""
Ideally the output JSON data files are in the same format as the Timur
consignments, which means we can re-use UI components.
This format converts a data frame / vector from


        col1    col2    col3
row1     1        3      5.1
row2     0       -8      87
row3     99      1       12.3

to an object like (note that multiple files should be
squashed together -- only the first matrix input
is shown here):

{
    "variable_name": {
        "matrix": {
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
    },
    "other_variable": {
        "vector": [
            {"label": "abc", "value": 1},
            {"label": "xyz", "value": 2}]
    },
    "variable_too": {
        "vector": [
            {"label": "123", "value": {
                "vector": [
                    {"label": "123.1", "value": 0.123},
                    {"label": "123.2", "value": 1.123}
                ]
            }},
            {"label": "124", "value": {
                "vector": [
                    {"label": "124.1", "value": 0.124},
                    {"label": "124.2", "value": 1.124}
                ]
            }}
        ]
    }
}

Don't know how nested vectors would look in a TSV, though...
"""

import argparse
import csv
import json
import tempfile
from pathlib import Path


class Matrix:
    def __init__(self, filename):
        pass

    @staticmethod
    def json():
        pass


class Vector:
    def __init__(self, filename):
        pass

    @staticmethod
    def json():
        pass


def isMatrix(csv_reader):
    """ Stupid test -- if row[0], col[0] is empty, is a matrix?
        otherwise is a label for a vector? """
    for row in csv_reader:
        if "" == row[0]:
            return False
        return True


def convert_tsv_to_json(inputs, output):
    with open(output, "w") as output_file:
        output_data = {}
        for input_filename in inputs:
            with open(input_filename, "r") as input_file:
                basename = Path(input_filename).stem
                reader = csv.reader(input_file, delimiter="\t")

                if isMatrix(reader):
                    data_class = Matrix(input_filename)
                else:
                    data_class = Vector(input_filename)
                file_output_data = {
                    "rows": []
                }
                headers = False
                row_names = []
                for row in reader:
                    row_data = row[1:]
                    if not headers:
                        headers = True
                        file_output_data["col_names"] = row_data
                        file_output_data["num_cols"] = len(row_data)
                        continue
                    row_names.append(row[0])

                    numeric_row = []
                    for datum in row_data:
                        try:
                            numeric_row.append(float(datum))
                        except:
                            numeric_row.append(datum)

                    file_output_data["rows"].append(numeric_row)
                file_output_data["num_rows"] = len(row_names)
                file_output_data["row_names"] = row_names

                output_data[basename] = file_output_data.copy()
        json.dump(output_data, output_file, indent=2)
        return output_data


def test_convert():
    results = convert_tsv_to_json(
        ["fixtures/data.tsv"], tempfile.NamedTemporaryFile().name)
    assert (results == {
        "data": {
            "matrix": {
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
            }
        }
    })


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert TSV to Consignment JSON format')
    parser.add_argument('--input', '-i', type=str,
                        default='data.tsv',
                        help='Input TSV file(s)',
                        nargs='*',
                        dest='inputs')
    parser.add_argument('--output', '-o', type=str,
                        default='data.json',
                        help='Output JSON file',
                        dest='output')

    arguments = parser.parse_args()
    convert_tsv_to_json(arguments.inputs, arguments.output)
