library(dataflow)

cluster_meta <- input_str('cluster_meta')
discrete_metadata_summary <- input_json('discrete_metadata_summary')

values <- as.list(discrete_metadata_summary[[cluster_meta]])
names(values) <- values

output_json(values, 'blank_annots.json')