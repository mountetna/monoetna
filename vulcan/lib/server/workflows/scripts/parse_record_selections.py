from archimedes.functions.dataflow import output_json, input_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import connect, question
from archimedes.functions.list import unique, flatten
from archimedes.functions.environment import project_name

pdat = input_json("project_data")[project_name]
selection_options = list(pdat['selection_options'].values())

magma = connect()

select1 = input_json('select1')
select2 = input_json('select2')
select3 = input_json('select3')

# Experiment and Tissue (AND logic)
filters = []
if len(select1) > 0:
    filters.append(
        buildTargetPath( selection_options[0], pdat ) + ['::in', select1]
    )
if len(select2) > 0:
    filters.append(
        buildTargetPath( selection_options[1], pdat ) + ['::in', select2]
    )
if len(select3) > 0:
    filters.append(
        buildTargetPath( selection_options[2], pdat ) + ['::in', select3]
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
