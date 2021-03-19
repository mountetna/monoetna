from archimedes.functions.dataflow import input_var, output_json

selected_type = input_var('type')

# These should really be fetched from Magma?
# The values should match what is used in user_filter_options.py
if selected_type == 'Experiment + Tissue':
    options = ['High-fat', 'Low-fat', 'Blood', 'Tissue']
elif selected_type == 'Pools':
    options = ['Pool1', 'Pool2']
else:
    options = [
        'XCRS1-MM170-SCPYMT3PTL1-SCG1',
        'XCRS1-MM170-SCPYMT3PTM1-SCG1'
    ]

output_json(options, 'options')
