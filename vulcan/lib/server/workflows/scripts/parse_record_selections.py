from archimedes.functions.dataflow import output_json, input_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.list import unique, flatten
from archimedes.functions.environment import project_name

pdat = input_json("project_data")[project_name]
selection_atts = pdat['selection_options']

# because a is the input name from the previous CWL input step
selected = input_json('selected_options')["a"]

magma = connect()

# Create filters for all the 'select-bys' that the value for this attribute
# must be among the options selected in the previous step.
filters = []
for selected_model in list(selected.keys()):
    if len(selected[selected_model]) > 0:
        filters.append(
            buildTargetPath( selection_atts[selected_model], pdat ) + ['::in', selected[selected_model]]
        )

seq_target = parseModelAttr(pdat['seq_h5_counts_data'])

tube_records = unique(question(
    magma,
    [
    seq_target['model'],
        [ '::has', seq_target['attribute']],
        *filters,
        '::all', '::identifier'
    ]
))

if len(tube_records) < 1:
    raise RuntimeError('No records with data meet the selected criteria.')

output_json(tube_records, 'tube_recs')
