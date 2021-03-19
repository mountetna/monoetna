from archimedes.functions.dataflow import output_path, input_path, output_json

mock_experiments = [
    'High-fat',
    'Low-fat'
]
output_json(mock_experiments, 'experiments')

mock_tissues = [
    'Blood',
    'Tumor'
]
output_json(mock_tissues, 'tissues')

mock_pools = [
    'Pool1',
    'Pool2'
]
output_json(mock_pools, 'pools')


mock_records = [
    'XCRS1-MM170-SCPYMT3PTL1-SCG1',
    'XCRS1-MM170-SCPYMT3PTM1-SCG1'
]
output_json(mock_records, 'records')
