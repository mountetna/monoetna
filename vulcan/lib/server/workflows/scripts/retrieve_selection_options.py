from archimedes.functions.dataflow import input_var, input_json, output_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique, flatten

pdat = input_json("project_data")
selection_options = pdat['selection_options']

seq_target = parseModelAttr(pdat['seq_h5_counts_data'])
q_start = [
    seq_target['model'],
    [ '::has', seq_target['attribute']],
    '::all'
]

magma = connect()

options = dict([
    [ 
        key,
        unique(flatten(question(
            magma,
            q_start + buildTargetPath( selection_options[key], pdat )
        )))
    ] for key in list(selection_options.keys()) ])

output_json(options, 'selection_options')
