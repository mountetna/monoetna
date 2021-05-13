from archimedes.functions.dataflow import output_json

project_data = {
	'xcrs1': {
		'seq_h5_counts_data': 'sc_seq#raw_counts_h5',  
		'seq_to_model_paths': {
			# Path should be what could go between '<model_of_above>' and '<attribute_name>' to magma/query any items from this model!
			'experiment': ['biospecimen_group', 'experiment'],
			'biospecimen_group': ['biospecimen_group'],
			'biospecimen': ['biospecimen_group', 'biospecimen'],
			'subject': ['biospecimen_group', 'biospecimen', 'subject'],
			'sc_seq_pool': ['sc_seq_pool'],
			'sc_seq': []
		},
		'color_options': {
			# Cluster, Tube, and Gene are standard and not needed here!
			# Format = <Label for the color-by drop down>: '<model>#<attribute>'
			'Experiment': 'experiment#alias',
			'Tissue': 'biospecimen_group#biospecimen_type',
			'Pool': 'sc_seq_pool#::identifier',
			'Biospecimen Group': 'biospecimen_group#::identifier',
		},
		'selection_options': {
			# Format = <Label of this selection item>: '<model>#<attribute>'
			'Experiment': 'experiment#alias',
			'Tissue': 'biospecimen_group#biospecimen_type',
			'Fraction': 'sc_seq#cell_fraction'
		}
	}
}

output_json(project_data, "project_data")