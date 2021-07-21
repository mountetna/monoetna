from archimedes.functions.dataflow import output_json

project_data = {
	'xcrs1': {
		'seq_h5_counts_data': 'sc_seq#raw_counts_h5',  
		'seq_to_model_paths': {
			# Path should be what could go between '<model_of_above>' and '<attribute_name_of_below>' to magma/query any items from that model!
			#   Great example: HuMu/xcrs1 'subject' where the path is to go...
			#   from 'sc_seq' to it's parent, 'biospecimen_group',
			# 		= ['biospecimen_group']
			#   then follow the link to 'biospecimen' (1:many; but note this link type-info is not actually needed!),
			# 		= ['biospecimen_group', 'biospecimen']
			#   then from there to travel up to that model's parent, 'subject'.
			# 		= ['biospecimen_group', 'biospecimen', 'subject']
			'experiment': ['biospecimen_group', 'experiment'],
			'biospecimen_group': ['biospecimen_group'],
			'biospecimen': ['biospecimen_group', 'biospecimen'],
			'subject': ['biospecimen_group', 'biospecimen', 'subject'],
			'sc_seq_pool': ['sc_seq_pool'],
			'sc_seq': []
		},
		'color_options': {
			# Format = <Label for the color-by drop down>: '<model>#<attribute>'
			#   Cluster, Tube, and Gene are standard options that do not need to be added here!
			'Experiment': 'experiment#alias',
			'Tissue': 'biospecimen_group#biospecimen_type',
			'Fraction': 'sc_seq#cell_fraction',
			'Pool': 'sc_seq_pool#::identifier',
			'Biospecimen Group': 'biospecimen_group#::identifier',
		},
		'selection_options': {
			# Format = <Label of this selection item>: '<model>#<attribute>'
			'Experiment': 'experiment#alias',
			'Tissue': 'biospecimen_group#biospecimen_type',
			'Fraction': 'sc_seq#cell_fraction'
		}
	},
	'ipi': {
		'seq_h5_counts_data': 'sc_rna_seq#raw_counts_h5',
		'seq_to_model_paths': {
			# Path should be what could go between '<model_of_above>' and '<attribute_name_of_below>' to magma/query any items from that model!
			#   Great example: HuMu/xcrs1 'subject' where the path is to go...
			#   from 'sc_seq' to it's parent, 'biospecimen_group',
			# 		= ['biospecimen_group']
			#   then follow the link to 'biospecimen' (1:many; but note this link type-info is not actually needed!),
			# 		= ['biospecimen_group', 'biospecimen']
			#   then from there to travel up to that model's parent, 'subject'.
			# 		= ['biospecimen_group', 'biospecimen', 'subject']
			'experiment': ['sample', 'patient', 'experiment'],
			'patient': ['sample', 'patient'],
			'sample': ['sample'],
			'sc_rna_seq': []
		},
		'color_options': {
			# Format = <Label for the color-by drop down>: '<model>#<attribute>'
			#   Cluster, Tube, and Gene are standard options that do not need to be added here!
			'Indication': 'experiment#name',
			'Tissue': 'sample#tissue_type',
			'Compartment': 'sc_rna_seq#biospecimen',
			'Chemistry': 'sc_rna_seq#chemistry',
			'Frozen tissue': 'patient#ffpe_frozen'
			# age, sex at birth, age at diag, bmi, smoker status, alcohol use, race, ethnicity, *past medical history (up to 10 per), time on ice
		},
		'selection_options': {
			# Format = <Label of this selection item>: '<model>#<attribute>'
			'Indication': 'experiment#name',
			'Tissue': 'sample#tissue_type',
			'Compartment': 'sc_rna_seq#biospeciman',
			'Chemistry': 'sc_rna_seq#chemistry'
		}
	}
}

output_json(project_data, "project_data")