from archimedes.functions.dataflow import output_json, input_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.list import unique, flatten

pdat = input_json("project_data")
selection_atts = pdat['selection_options']

selected = input_json('selected_options')

magma = connect()

# Create filters for all the 'select-bys' that the value for this attribute
# must be among the options selected in the previous step.
filters = []
for target in list(selection_atts.keys()):
    if len(selected[target]) > 0:
        filters.append(
            buildTargetPath( selection_atts[target], pdat ) + ['::in', selected[target]]
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
