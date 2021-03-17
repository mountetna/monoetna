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
