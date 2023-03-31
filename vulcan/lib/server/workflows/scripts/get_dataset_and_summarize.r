suppressPackageStartupMessages({
    library(Seurat)
    library(magmaR)
    library(dataflow)
    library(environment)
    library(dittoSeq)
})

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
    lapply(names, USE.NAMES = TRUE, FUN = function(x) {
        as.character(seq_len(ncol(
            embeds_getter(object, x)
        )))
    })
} 
summarize_plotting_options <- function(object, record) {
    # Outputs:
    # - plotting_options = helps the adaptor function convert user selections into workable inputs to the functions

    metadata_names <- getMetas(object)

    # Cell Metadata & Dimensionality Reduction Embeddings
    plotting_options = list(
        'Cell_Metadata' = getMetas(object)
        'Reductions' = getReductionsDims(object)
    )
    # Features
    if (!is.na(record["modalities"])) {
        modalities <- unlist(strsplit(record["modalities"], ","))
        plotting_options$Cell_Features <- list()
        if ("gex" %in% modalities) {
            plotting_options$Cell_Features$RNA <- getGenes(object, assay = "RNA")
        }
        if ("cite" %in% modalities) {
            plotting_options$Cell_Features$ADT <- getGenes(object, assay = "ADT")
        }
    }

    ### Recommendations = bits marked in the record as putative
    clust_and_annot_vals <- record[c("metadata_clustering", "metadata_annots_fine", "metadata_annots_broad")]
    if (any(!is.na(clust_and_annot_vals))) {
        names(clust_and_annot_vals) <- c("clustering", "annotations_fine", "annotations_broad")
        plotting_options["Recommended_Metadata"] <- as.list(clust_and_annot_vals[!is.na(clust_and_annot_vals)])
    }
    if (!is.na(record["umap_name"])) {
        plotting_options["Recommended_Reduction"] <- record["umap_name"]
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

            output[[i]] <- as.list(this_out)
        }
    }

    if (!is.null(recommendations)) {
        if (!is.null(recommendations$clustering)) output$`_Clustering_` <- outputs[[recommendations$clustering]]
        if (!is.null(recommendations$annotations_fine)) output$`_Annotations_Fine_` <- outputs[[recommendations$annotations_fine]]
        if (!is.null(recommendations$annotations_broad)) output$`_Annotations_Broad_` <- outputs[[recommendations$annotations_broad]]
    }

    output
}

name_wrap <- function(vector) {
    new_vals <- rep(NULL, length.out = length(vector))
    names(new_vals) <- vector
    as.list(new_vals)
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
}

magma = magmaRset(
    token = environment::token(),
    url = environment::magma_host()
)

dataset_record_raw <- retrieve(
    magma,
    project_name(),
    "sc_seq_dataset",
    recordNames = input_str(dataset_name),
    attributeNames = "all"
)
dataset_record <- as.vector(dataset_record_raw)
names(dataset_record) <- colnames(dataset_record_raw)

# Get the data
dataset_url <- dataset_record["object"]
# Download directly to output location!
dataset_file <- HttpClient$new(url = dataset_url)$get(disk = output_path("scdata")) # (Outputs the filepath)
scdata = readRDS(dataset_file)

plotting_options <- summarize_plotting_options(scdata, dataset_record)
discrete_metadata_summary <- summarize_discrete_metadata(scdata)

reduction_opts <- plotting_options[['Reductions']] # [[]] returns the internal value
if (!is.null(plotting_options["Recommended_Reduction"])) {
    reduction_opts$`_Recommended_` <- reduction_opts[[plotting_options["Recommended_Reduction"]]]
}
continuous_opts <- plotting_options['Cell_Features'] # [] subsets the top-level list to the called elements, keeping the top-level
continuous_opts$Cell_Metadata <- setdiff(plotting_options[['Cell_Metadata']], names(discrete_metadata_summary))
discrete_opts <- list('Cell_Metadata' = names(discrete_metadata_summary))
all_opts <- plotting_options[c('Cell_Features', 'Cell_Metadata')]

# Outputs
# dataset downloaded to final location already
output_json(plotting_options, 'plotting_options')
output_json(discrete_metadata_summary, 'discrete_metadata_summary')
output_json(toNestedOptionsSet(all_opts), 'all_opts')
output_json(toNestedOptionsSet(continuous_opts), 'continuous_opts') # Needs to be nestedOptionSet
output_json(toNestedOptionsSet(discrete_opts), 'discrete_opts')
output_json(reduction_opts, 'reduction_opts')
