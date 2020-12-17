###### magmaR Demo - Nov 19, 2020 ######

### Installation
# remotes::install_github("mountetna/etna", subdir = "packages/magmaR")
# Any dependencies will be installed, automatically, based on contents of
# DESCRIPTION file `Imports:` field.

library(magmaR)
?retrieve

# Operates via 'RCurl', with json handling via 'jsonlite', and string-tsv -> dataframe conversion via 'readr' 

### TOKEN utilization:
# Obtained in first call via prompt to user
ids <- retrieveIds(
    projectName = "ipi",
    modelName = "rna_seq")
# Saved as a hidden object, .MAGMAR_TOKEN.
.MAGMAR_TOKEN

# Can also be given explicitly:
ids_patient <- retrieveIds(
    projectName = "ipi",
    modelName = "patient",
    token = "<give_token_here>")

# Error message when magma sends back that user is unauthorized:
#   You are unauthorized. If you think this is a mistake, run `rm(.MAGMAR_TOKEN)` or update your 'token' input, then retry. 

### Dev / Staging / Production selection:
# 'url.base' input can be given to any function.
ids_patient <- retrieveIds(
    projectName = "ipi",
    modelName = "patient")

# When not given, what's used is based on a '.MAGMAR_URL' variable which is created upon first function call.
# Default is "https://magma.ucsf.edu"
.MAGMAR_URL

# To update this default value for when 'url.base' is not given, run code like this:
.MAGMAR_URL <- "http://magma.development.local" # = development

### Main functions:

## "basic" /query
query(
    projectName = "ipi",
    queryTerms = 
        list('rna_seq',
             '::all',
             'sample',
             '::identifier') ### Note from Saurabh: /query can also pull multiple attributes, so making a pretty df may be more work than currently implemented.
    # , format = "df" ### Note from Saurabh: format may have more than just column names, so need more info to assess consistent df output viability.
    )
# ^^^: output = list by default, but can be dataframe as well

## "basic" /retrieve
df <- retrieve(
    projectName = "ipi",
    modelName = "patient",
    recordNames = "all",
    attributeNames = "all",
    filter = "ipi_number~GYN", # filter works =)
    pageSize = 25,
    page = 1
    )
# ^^^: 'format = "tsv"', output = dataframe

json <- retrieveJSON(
    projectName = "ipi",
    modelName = "patient",
    recordNames = "all",
    attributeNames = "all")
# ^^^: 'format = "json"', output = list

## "Specialized" /retrieve

mat <- retrieveMatrix(
    projectName = "ipi",
    modelName = "rna_seq",
    recordNames = ids[grepl("live2",ids)],
    attributeNames = "gene_tpm")
# ^^^: Under the hood, .matrix_retrieval_chunk() grabs data as a json for 10 records at a time.
# These are then converted to df columns & cbind'd, and rownames are grabbed from the template.

meta <- retrieveMetadata(
    projectName = "ipi",
    meta_modelName = "population",
    meta_attributeNames = "all",
    target_modelName = "rna_seq",
    target_recordNames = ids[grepl("live2",ids)])
# ^^^: internals which I think give a good idea of how it's done:
# .trace_model_to_proj (out = vector of modelNames)
# .obtain_linkage_paths (out = 1 .trace, or 2 .traces trimmed to first intersect, where '.trace'=output of the fxn above)
# .map_identifiers_by_path (out = data frame of recordNames with columns for each model of a path; run separately for target data & for metadata if models are on separate branches)
# .expand_metadata_to_have_1row_per_id (out = cbind'ed retrieve() data slices; run ONLY when target and meta models are on separate branches and so can have 1:many id mapping)

### Helper functions:

# To determine `modelName`, recordNames` or `attributeNames` options:
retrieveModels(
    projectName = "ipi")

retrieveIds(
    projectName = "ipi",
    modelName = "patient")

retrieveAttributes(
    projectName = "ipi",
    modelName = "patient")

# To retrieve the project template:
temp <- retrieveTemplate(
    projectName = "ipi")

#### Plans:
# 1) Build out unit testing. # webmockr / vcr
# 2) Create vignette.
# 3) Finalize inputs / names of current functions.
#
# 4) Release to DS team?
#
# DONE 5) Better dev / staging integration to remove need for `url.base` explicitness?
# 6) Add a 'group.by' option to 'retrieveMetadata()'.
# 7) Implement 'connected.only' in `retrieve()`/`query()` once ready.
