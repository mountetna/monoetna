from archimedes.functions.dataflow import output_json, input_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.list import unique, flatten
from archimedes.functions.environment import project_name

pdat = input_json("project_data")[project_name]
selection_atts = pdat['selection_options']

selected = input_json('selected_options')

magma = connect()

# Experiment and Tissue (AND logic)
filters = []
for set in list(selected.keys()):
    if len(selected[set]) > 0:
        filters.append(
            buildTargetPath( selection_atts[set], pdat ) + ['::in', selected[set]]
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

output_json(tube_records, 'tube_recs')
