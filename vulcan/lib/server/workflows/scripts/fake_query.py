from archimedes.functions.dataflow import output_path, input_path, output_json

experiments = [
    'all',
    'XCRS1-SCPYMT3'
]

tissues = [
    'all',
    'Tumor',
    'Blood',
    'LN'
]

pools = [
    ''
]

recs = [
    ''
    'XCRS1-MM170-SCPYMT3PTL1-SCG1',
    'XCRS1-MM170-SCPYMT3PTM1-SCG1'
]

output_json(experiments, 'experiments')
output_json(tissues, 'tissues')
output_json(pools, 'pools')
output_json(recs, 'records')
