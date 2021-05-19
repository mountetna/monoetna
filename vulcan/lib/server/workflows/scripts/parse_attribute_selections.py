from archimedes.functions.dataflow import output_json, input_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.list import unique, flatten
from archimedes.functions.environment import project_name

scdata = sc.read(input_path('scdata'))
pdat = input_json("project_data")[project_name]

# because a is the input name from the previous CWL input step
selected = input_json('selected_options')["a"]

# Build actual options
magma = connect()

seq_target = parseModelAttr(pdat['seq_h5_counts_data'])
q_start = [
    seq_target['model'],
    [ '::identifier', '::in', scdata.obs[ 'Record_ID' ]],
    '::all'
]
s
election_values = dict([
    att,
    unique(flatten(question(
        magma,
        [
            q_start + buildTargetPath( selection_atts[target], pdat )
        ]
    )))
] for att in selected)

output_json(selection_values, 'selection_values')
