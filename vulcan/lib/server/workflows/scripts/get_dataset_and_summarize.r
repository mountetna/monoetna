library(Seurat)
library(magmaR)
library(dataflow)
library(dittoSeq)

# A Re-implementation of dittoSeq's getReductions because we will also need to know how many components each reductions has
getReductionsDims <- function(object) {
    # Return:
    # named list object where:
    #   names are the reduced dims that exist within the object
    #   values are 1:#dims in each reduced dim.

    if (is(object,"SingleCellExperiment")) {
        name_getter <- SingleCellExperiment::reducedDimNames
        embeds_getter <-SingleCellExperiment::reducedDim
    }
    # Seurat modern
    if (is(object,"Seurat")) {
        name_getter <- embeds_getter <- Seurat::Reductions
    }
    # # Seurat v2 (DL won't need soon. Leaving as task for when moving to dittoSeq)
    # if (is(object,"seurat")) {
    #     return(names(object@dr))
    # }
    names <- name_getter(object)
    out <- lapply(names, FUN = function(x) {
        as.character(seq_len(ncol(
            embeds_getter(object, x)
        )))
    })
    names(out) <- names
    out
} 
summarize_plotting_options <- function(object, record) {
    # Outputs:
    # - plotting_options = helps the adaptor function convert user selections into workable inputs to the functions

    metadata_names <- getMetas(object)

    # Cell Metadata & Dimensionality Reduction Embeddings
    plotting_options = list(
        'Cell_Metadata' = getMetas(object),
        'Reductions' = getReductionsDims(object)
    )
    # Features
    if (!is.na(record[["modalities"]])) {
        modalities <- unlist(strsplit(record[["modalities"]], ","))
        plotting_options$Cell_Features <- list()
        if ("gex" %in% modalities) {
            plotting_options$Cell_Features$RNA <- getGenes(object, assay = "RNA")
        }
        if ("cite" %in% modalities) {
            plotting_options$Cell_Features$ADT <- getGenes(object, assay = "ADT")
        }
    }

    ### Recommendations = bits marked in the record as putative
    meta_recs <- c()
    for (rec in c("metadata_clustering", "metadata_annots_fine", "metadata_annots_broad")) {
        if (!is.null(record[[rec]])) {
            new <- record[[rec]]
            names(new) <- paste0(gsub("^metadata_", "_Recommended_", rec), "_")
            meta_recs <- c(meta_recs, new)
        }
    }
    if (length(meta_recs) > 0) {
        plotting_options[["Recommended_Metadata"]] <- as.list(meta_recs)
    }
    if (!is.na(record[["umap_name"]])) {
        plotting_options[["Recommended_Reduction"]] <- record[["umap_name"]]
    }

    plotting_options
}

summarize_discrete_metadata <- function(
    object,
    recommendations = NULL # If given, should be `summarize_plotting_options()$Recommended_Metadata`
    ) {
    
    metadata_values <- as.data.frame(getMetas(scdata, names.only = FALSE))
    
    output <- list()
    for (i in getMetas(scdata)) {
        this_data <- metadata_values[,i, drop = TRUE]
        if (!is.numeric(this_data)) {
            # Remove any empty levels, but keep order if already a factor
            used <- levels(as.factor(as.character(this_data)))
            levs <- levels(as.factor(this_data))
            this_out <- levs[levs %in% used]

            # 'key' as numbers to match with the pandas output format already being used
            names(this_out) <- seq_along(this_out)

            output[[i]] <- this_out
        }
    }

    if (!is.null(recommendations)) {
        if (!is.null(recommendations$`_Recommended_clustering_`)) output$`_Recommended_clustering_` <- output[[recommendations$`_Recommended_clustering_`]]
        if (!is.null(recommendations$`_Recommended_annots_fine_`)) output$`_Recommended_annots_fine_` <- output[[recommendations$`_Recommended_annots_fine_`]]
        if (!is.null(recommendations$`_Recommended_annots_broad_`)) output$`_Recommended_annots_broad_` <- output[[recommendations$`_Recommended_annots_broad_`]]
    }

    output
}

name_wrap <- function(vector) {
    new_vals <- lapply(seq_along(vector), function(x) {return(NULL)})
    names(new_vals) <- vector
    new_vals
}
toNestedOptionsSet <- function(named_list) {
    # Finds all leaves in the list that are not (NULL | named list) and converts them to a named list with value of null
    # = Formatting needed for nestedDropdown UIs
    for (i in seq_along(named_list)) {
        if (!is.null(named_list[[i]])) {
            named_list[[i]] <- 
                if (!is.list(named_list[[i]])) {
                    name_wrap(named_list[[i]])
                } else {
                    toNestedOptionsSet(named_list[[i]])
                }
        }
    }
    named_list
}

magma_opts <- list(
    followlocation = FALSE)
if (grepl("development", magma_host())) {
    magma_opts$ssl_verifyhost <- FALSE
    magma_opts$ssl_verifypeer <- FALSE
}
magma = magmaRset(
    token = token(),
    url = magma_host(),
    opts = magma_opts
)

dataset_record_raw <- retrieve(
    magma,
    project_name(),
    "sc_seq_dataset",
    recordNames = input_str('dataset_name'),
    attributeNames = "all"
)
dataset_record <- as.vector(dataset_record_raw)
names(dataset_record) <- colnames(dataset_record_raw)

# Get the data
dataset_url <- dataset_record[["object"]]
dataset_path <- output_path("scdata")
# Download directly to output location!
x <- crul::HttpClient$new(url = dataset_url)$get(disk = dataset_path) # Outputs status
scdata <- readRDS(dataset_path)

plotting_options <- summarize_plotting_options(scdata, dataset_record)

discrete_metadata_summary <- summarize_discrete_metadata(scdata, plotting_options$Recommended_Metadata)

reduction_opts <- plotting_options[['Reductions']] # [[]] returns the internal value
if (!is.null(plotting_options[["Recommended_Reduction"]])) {
    reduction_opts$`_Recommended_` <- reduction_opts[[plotting_options[["Recommended_Reduction"]]]]
}

continuous_opts <- plotting_options['Cell_Features'] # [] subsets the top-level list to the called elements, keeping the top-level
continuous_opts$Cell_Metadata <- setdiff(plotting_options[['Cell_Metadata']], names(discrete_metadata_summary))

discrete_opts <- names(discrete_metadata_summary)

all_opts <- plotting_options[c('Cell_Features', 'Cell_Metadata')]

# Outputs
# dataset downloaded to final location already
output_json(plotting_options, 'plotting_options')
output_json(discrete_metadata_summary, 'discrete_metadata_summary')
output_json(toNestedOptionsSet(all_opts), 'all_opts')
output_json(toNestedOptionsSet(continuous_opts), 'continuous_opts') # Needs to be nestedOptionSet
output_json(discrete_opts, 'discrete_opts')
output_json(reduction_opts, 'reduction_opts')
