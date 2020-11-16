library(magmaR)
library(jsonlite)

# Attribute Names
ipi_rnaseq_atts <- retrieveAttributeNames("ipi", "rna_seq")

# List / (ideally dataframe, but they aren't rectangular) retrieval
rnaseq_data <- retrieve("ipi", "rna_seq", attributeNames = ipi_rnaseq_atts[-1*c(15,16)])
rnaseq_data <- retrieve("ipi", "rna_seq", attributeNames = "all")
rnaseq_data <- retrieve("ipi", "rna_seq", attributeNames = ipi_rnaseq_atts[1:2])

# Identifiers
ipi_rnaseq_ids <- retrieve("ipi", "rna_seq", attributeNames = "identifier")

# Retrieve a single one
ipi_rnaseq_rec1 <- retrieveJSON("ipi", "rna_seq", recordNames = "IPIADR002.T1.rna.live", attributeNames = "gene_tpm")
####### ERROR

# gene counts
ipi_rnaseq_rec1$models$rna_seq$documents$IPIADR002.T1.rna.live$gene_tpm
# gene names
ipi_rnaseq_rec1$models$rna_seq$template$attributes$gene_tpm$options

# Retrieve multiple
ipi_rnaseq_10recs <- retrieveJSON("ipi", "rna_seq", recordNames = ipi_rnaseq_ids[40:49], attributeNames = "gene_counts")

# Retrieve all gene_counts
# DO NOT RUN THIS. retrieveMatrix was made for this exact purpose... this json is to big to handle.
# Reminder to test though & maybe add some sort of exception...
ipi_rnaseq_all_counts <- retrieveJSON("ipi", "rna_seq", recordNames = "all", attributeNames = "gene_counts")

# Retrieve 1 sample's gene_tpm
ipi_rnaseq_samp1_counts <- retrieveJSON(
    "ipi", "rna_seq", recordNames = "IPIADR002.T1.rna.live", attributeNames = "gene_counts")

# retrieve page 1
(ipi_rnaseq_page1 <- retrieveJSON(
    "ipi", "rna_seq", recordNames = "IPIADR002.T1.rna.live", attributeNames = "gene_counts",
    page = 1))
# > ERROR

# Page retrieval
(ipi_rnaseq_page1 <- retrieveJSON(
    "ipi", "rna_seq", recordNames = "all", attributeNames = "compartment",
    page = 1))
# > ERROR

(ipi_rnaseq_page2 <- retrieveJSON(
    "ipi", "rna_seq", recordNames = "all", attributeNames = "compartment",
    page = 2))
# > ERROR

# Just get 1 full rnaseq record
IPIADR002.T1.rna.live <- retrieveJSON(
    "ipi", "rna_seq", recordNames = "IPIADR002.T1.rna.live", attributeNames = "gene_counts")

recs2 <- retrieveJSON(
    "ipi", "rna_seq",
    recordNames = c("IPIADR002.T1.rna.live", "IPIADR002.T1.rna.stroma"),
    attributeNames = "gene_counts")

id <- 7
retrieveJSON(
    "ipi", "rna_seq", recordNames = ipi_rnaseq_ids[id], attributeNames = "gene_counts")$models$rna_seq$documents[[ipi_rnaseq_ids[id]]]$gene_counts[1:10]

# If BUILT
magmaR:::.matrix_retrieval_chunk(
    "ipi", "rna_seq", recordNames = ipi_rnaseq_ids[1:10], attributeNames = "gene_counts")
# If SOURCED
.matrix_retrieval_chunk(
    "ipi", "rna_seq", recordNames = ipi_rnaseq_ids[1:10], attributeNames = "gene_counts")

ipi_rnaseq_ids <- retrieve("ipi", "rna_seq", attributeNames = "identifier")
rnaseq_1to100 <- retrieveMatrix(
    "ipi", "rna_seq", recordNames = ipi_rnaseq_ids[1:100], attributeNames = "gene_counts")

rnaseq_1to100_json <- toJSON(rnaseq_1to100)
# NOTE = 19.6 mb

rnaseq_live <- retrieveMatrix(
    "ipi", "rna_seq", recordNames = ipi_rnaseq_ids[grep("live",ipi_rnaseq_ids)], attributeNames = "gene_counts")
dim(rnaseq_live)
rnaseq_live[1:10,1:10]

rnaseq_live_json <- toJSON(rnaseq_live)
# NOTE, takes ~5sec and 133.6 mb!!

temp <- retrieveTemplate("ipi")


#############
### query ###
#############

samples_of_rna_seq <- query(
    "ipi",
    queryTerms = 
        list('rna_seq',
             list('::has', 'gene_tpm'),
             '::all',
             'sample',
             '::identifier'
        ))


##################
### id mapping ###
##################

temp <- retrieveTemplate("ipi")
.trace_model_to_proj("rna_seq", temp)
# Gives: "rna_seq"    "sample"     "patient"    "experiment" "project" 

(paths <- .obtain_linkage_paths("rna_seq", "treatment", temp))
# Gives:
# $target_path
# [1] "rna_seq" "sample"  "patient"
# 
# $meta_path
# [1] "treatment" "patient"

(target_id_map <- .map_identifiers_by_path(
    path = paths$target_path,
    projectName = "ipi"))
# df with 3507 rows and 3 cols
(meta_id_map <- .map_identifiers_by_path(
    path = paths$meta_path,
    projectName = "ipi"))
# df with 771 rows and 2 cols

ids <- retrieveIds("ipi", "rna_seq")
meta <- retrieveMetadata(
    projectName = "ipi",
    target_modelName = "rna_seq",
    target_recordNames = ids[grepl("live", ids)],
    meta_modelName = "treatment",
    meta_attributeNames = "all")
# GOOD
    
meta <- retrieveMetadata(
    projectName = "ipi",
    target_modelName = "rna_seq",
    target_recordNames = ids[grepl("live", ids)],
    meta_modelName = "sample",
    meta_attributeNames = "all")
# GOOD

meta <- retrieveMetadata(
    projectName = "ipi",
    target_modelName = "rna_seq",
    target_recordNames = "all",
    meta_modelName = "sample",
    meta_attributeNames = "all")
# GOOD

should_be_df <- retrieve(
    projectName = "ipi",
    modelName = "sample",
    recordNames = "all",
    attributeNames = "all")

should_be_df <- retrieve(
    projectName = "ipi",
    modelName = "rna_seq",
    recordNames = ids[grepl("live3",ids)],
    attributeNames = "all")

should_be_matrix <- retrieveMatrix(
    projectName = "ipi",
    modelName = "rna_seq",
    recordNames = ids[grepl("live2",ids)],
    attributeNames = "gene_tpm")

json <- toJSON(as.data.frame(should_be_matrix))
df <- fromJSON(json)
rownames(df)
matrix <- as.matrix(df)


