library(dataflow)

plotting_options <- input_json('plotting_options')
discrete_metadata_summary <- input_json('discrete_metadata_summary')

opts <- names(discrete_metadata_summary)[names(discrete_metadata_summary) %in% plotting_options$Cell_Metadata]
rec <- if ('_Recommended_clustering_' %in% names(plotting_options[["Recommended_Metadata"]])) {
    plotting_options$Recommended_Metadata$`_Recommended_clustering_`
} else {
    NULL
}

output_string(rec, 'rec')
output_json(opts, 'opts')
