from archimedes.functions.dataflow import output_path, input_path, output_json

mock_data = [
    'XCRS1-MM170-SCPYMT3PTL1-SCG1',
    'XCRS1-MM170-SCPYMT3PTM1-SCG1'
]
output_json(mock_data, 'tube_recs')
