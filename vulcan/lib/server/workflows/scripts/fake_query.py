from archimedes.functions.dataflow import output_path, input_path, json

with open(output_path('names'), 'w') as output_file:
    mock_data = [
        'XCRS1-MM170-SCPYMT3PTL1-SCG1',
        'XCRS1-MM170-SCPYMT3PTM1-SCG1'
    ]
    json.dump(mock_data, output_file, indent=2)
