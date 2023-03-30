suppressPackageStartupMessages({
    library(Seurat)
    library(magmaR)
    library(dataflow)
    library(environment)
    library(dittoSeq)
})

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
# Download directly to output location
dataset_file <- HttpClient$new(url = dataset_url)$get(disk = output_path("sc_data")) # (Outputs the filepath)
scdata = readRDS(dataset_file)

metadata_names <- getMetas(scdata)
metadata_values <- as.data.frame(getMetas(scdata, names.only = FALSE))

### Parse data in the object which might get used
# Cell Metadata & Dimensionality Reduction Embeddings
plotting_options = list(
    'Cell_Metadata' = metadata_names
    'Reductions' = getReductions(scdata)
)
# Features
if (!is.na(dataset_record["modalities"])) {
    plotting_options$Cell_Features <- list()
    if ("gex" %in% dataset_record["modalities"]) {
        plotting_options$Cell_Features$RNA <- getGenes(scdata, assay = "RNA")
    }
    if ("cite" %in% dataset_record["modalities"]) {
        plotting_options$Cell_Features$ADT <- getGenes(scdata, assay = "ADT")
    }
}

### Recommendations = bits marked in the record as putative
clust_and_annot_vals <- dataset_record[c("metadata_clustering", "metadata_annots_fine", "metadata_annots_broad")]
if (any(!is.na(clust_and_annot_vals))) {
    names(clust_and_annot_vals) <- c("clustering", "annotations_fine", "annotations_broad")
    plotting_options["Recommended_Metadata"] <- as.list(clust_and_annot_vals[!is.na(clust_and_annot_vals)])
}
if (!is.na(dataset_record["umap_name"])) {
    plotting_options["Recommended_Reduction"] <- dataset_record["umap_name"]
}

### Summarize metadata (b/c loading the entire data.frame in the next step, the users' browser, is not seem scalable.)
continuous_opts <- c()
discrete_opts <- c()
metadata_summary <- lapply(metadata_names, function(i) {
    this_data <- metadata_values[,i, drop = TRUE]
    if (is.numeric(this_data)) {
        continuous_opts <- c(continuous_opts, i)
        # (Really just a stub as these values will only be checked for discrete data)
        out <- c(min(this_data, na.rm = TRUE), max(this_data, na.rm = TRUE))
    } else {
        discrete_opts <- c(discrete_opts, i)
        # Remove any empty levels, but keep order if already a factor
        used <- levels(as.factor(as.character(this_data)))
        levs <- levels(as.factor(this_data))
        
        out <- levs[levs %in% used]
    }
    # 'key' as numbers to match with the pandas output format already being used
    names(out) <- seq_along(out)
    
    as.list(out)
})
names(metadata_summary) <- metadata_names

# Outputs
# dataset downloaded to final location already
output_json(plotting_options, 'plotting_options')
output_json(metadata_summary, 'metadata_summary')
output_json(continuous_opts, 'continuous_opts')
output_json(discrete_opts, 'discrete_opts')