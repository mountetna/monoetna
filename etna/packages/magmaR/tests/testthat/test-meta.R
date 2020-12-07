# This code tests the retrieveMetadata function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-meta.R")

test_that("retrieveMetadata, same branch", {
    vcr::use_cassette("metadata_same_branch", {
        atts <<- retrieveAttributes("example", "biospecimen")
        
        # metadata from distinct branch of project tree
        meta <-  retrieveMetadata(
            projectName = "example",
            target_modelName = "rna_seq",
            target_recordNames = "all",
            meta_modelName = "biospecimen",
            meta_attributeNames = "all")
    })
    
    expect_s3_class(meta, "data.frame")
    expect_gte(ncol(meta), length(atts))
    expect_equal(nrow(meta), 12)
    expect_lte(sum(!atts %in% colnames(meta)), 1)
})

# # Currently broken because of incomplete function logic:
# #   The biospecimen modl is the first shared model here, but the subject model is the first where
# #   Ids are actually shared between flow and rna_seq.
# test_that("retrieveMetadata, different branch", {
#     vcr::use_cassette("metadata_diff_branch", {
#         
#         # metadata from distinct branch of project tree
#         meta <- retrieveMetadata(
#             projectName = "example",
#             target_modelName = "rna_seq",
#             target_recordNames = "all",
#             meta_modelName = "flow",
#             meta_attributeNames = "all")
#     })
#     
#     expect_s3_class(meta, "data.frame")
#     expect_gt(ncol(meta), length(atts))
#     expect_equal(nrow(meta), 12)
#     expect_lte(sum(!atts %in% colnames(meta)), 1)
# })
