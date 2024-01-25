library(dataflow)

## This path is more complicated in code, but FAR FASTER than if we just added the metadata, saved the object, and then reran the summarization script.

name_wrap <- function(vector) {
    new_vals <- lapply(seq_along(vector), function(x) {return(NULL)})
    names(new_vals) <- vector
    new_vals
}

annots <- read.csv(input_path("annots.csv"), header = TRUE)
clus_meta <- names(annots)[1]
annot_metas <- grep("^annot", names(annots)[-1], value = TRUE)

plotting_options <- input_json('plotting_options')
discrete_metadata_summary <- input_json('discrete_metadata_summary')
all_opts <- input_json('all_opts')
discrete_opts <- input_json('discrete_opts')

# plotting_options & discrete_opts & all_opts
plotting_options$Cell_Metadata <- c(plotting_options$Cell_Metadata, annot_metas)
discrete_opts <- c(discrete_opts, annot_metas)
all_opts$Cell_Metadata <- c(all_opts$Cell_Metadata, name_wrap(annot_metas))

# discrete_metadata_summary
for (annot_meta in annot_metas) {
    this_meta_vals <- list(a = unique(annots[,annot_meta]))
    names(this_meta_vals) <- annot_meta
    discrete_metadata_summary[[annot_meta]] <- this_meta_vals
}

# Output
output_json(plotting_options, 'plotting_options')
output_json(discrete_metadata_summary, 'discrete_metadata_summary')
output_json(all_opts, 'all_opts')
output_json(discrete_opts, 'discrete_opts')
