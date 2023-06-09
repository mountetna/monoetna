from archimedes.functions.dataflow import input_var, input_json, output_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique, flatten

seq_target = parseModelAttr('sc_seq_dataset#name')
magma = connect()

options = unique(flatten(
    question(
        magma,
            [
                seq_target['model'],
                [ '::has', seq_target['attribute']],
                '::all'
            ]
        )
    ))

output_json(options, 'dataset_options')
