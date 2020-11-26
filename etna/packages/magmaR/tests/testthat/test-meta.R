# This code tests the retrieveMetadata function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-meta.R")

.GlobalEnv$.MAGMAR_TOKEN <- Sys.getenv("TOKEN")

test_that("retrieveMetadata, same branch", {
    vcr::use_cassette("metadata_same_branch", {
        ids <- retrieveIds("ipi", "rna_seq")
        atts <- retrieveAttributes("ipi", "sample")
        
        # metadata from distinct branch of project tree
        meta <-  retrieveMetadata(
            projectName = "ipi",
            target_modelName = "rna_seq",
            target_recordNames = ids[grep("live", ids)[1:100]],
            meta_modelName = "sample",
            meta_attributeNames = "all")
    })
    
    expect_s3_class(meta, "data.frame")
    expect_gte(ncol(meta), length(atts))
    expect_equal(nrow(meta), 100)
    expect_lte(sum(!atts %in% colnames(meta)), 1)
})

test_that("retrieveMetadata, different branch", {
    vcr::use_cassette("metadata_diff_branch", {
        ids <- retrieveIds("ipi", "rna_seq")
        atts <- retrieveAttributes("ipi", "demographic")
        
        # metadata from distinct branch of project tree
        meta <- retrieveMetadata(
            projectName = "ipi",
            target_modelName = "rna_seq",
            target_recordNames = ids[grep("live", ids)[1:100]],
            meta_modelName = "demographic",
            meta_attributeNames = "all")
    })
    
    expect_s3_class(meta, "data.frame")
    expect_gt(ncol(meta), length(atts))
    expect_equal(nrow(meta), 100)
    expect_lte(sum(!atts %in% colnames(meta)), 1)
})
