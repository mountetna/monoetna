from archimedes.functions.dataflow import output_path, input_path, json

with open(output_path('names'), 'w') as output_file:
    mock_data = [
        'XCRS1-MM431',
        'XCRS1-MM443',
        'XCRS1-MM204',
        'XCRS1-MM223'
    ]
    json.dump(mock_data, output_file, indent=2)
