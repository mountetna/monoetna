from archimedes.functions.dataflow import output_json, input_json, output_var

pdat = input_json("project_data")

no_batch='NO BATCH CORRECTION'
opts = [no_batch, 'Record_ID'] + list(pdat['color_options'].keys())

output_json(opts, 'batch_options')
output_var(no_batch, 'no_batch_string')