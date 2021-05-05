from archimedes.functions.dataflow import input_var, input_json, output_json, parseModelAttr, buildTargetPath
from archimedes.functions.magma import question, connect
from archimedes.functions.list import unique, flatten
from archimedes.functions.environment import project_name

pdat = input_json("project_data")[project_name]
selection_options = list(pdat['selection_options'].values())

seq_target = parseModelAttr(pdat['seq_h5_counts_data'])
q_start = [
    seq_target['model'],
    [ '::has', seq_target['attribute']],
    '::all'
]

magma = connect()

select1 = question(
    magma,
    q_start + buildTargetPath( selection_options[0], pdat )
)

select2 = question(
    magma,
    q_start + buildTargetPath( selection_options[1], pdat )
)

select3 = question(
    magma,
    q_start + buildTargetPath( selection_options[2], pdat )
)

output_json(unique(select1), 'select1')
output_json(unique(select2), 'select2')
output_json(unique(select3), 'select3')
