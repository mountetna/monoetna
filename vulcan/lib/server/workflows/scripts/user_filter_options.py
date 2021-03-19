from archimedes.functions.dataflow import output_json

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

color_options = {
    'Experiment': {
        'High-fat': None,
        'Low-fat': None
    },
    'Tissue': {
        'Blood': None,
        'Tissue': None
    },
    'Gene': {
        'Gene1': None,
        'Gene2': None,
        'Gene3': None
    }
}

output_json(experiments, 'experiments')
output_json(tissues, 'tissues')
output_json(color_options, 'color_options')
