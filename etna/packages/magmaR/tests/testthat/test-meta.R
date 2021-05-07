# This code tests the retrieveMetadata function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-meta.R")

targ <- magmaRset(
    token = TOKEN,
    url = URL)

test_that("retrieveMetadata, same branch", {
    vcr::use_cassette("metadata_same_branch", {
        atts <- retrieveAttributes(targ, "example", "biospecimen")
        
        # metadata from distinct branch of project tree
        meta <-  retrieveMetadata(
            target = targ,
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

# # Currently broken because there is no flow data loaded.
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
