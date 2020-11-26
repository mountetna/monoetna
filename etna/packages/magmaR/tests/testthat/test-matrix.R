# This code tests the retrieveMatrix function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-matrix.R")

test_that("retrieveMatrix", {
    vcr::use_cassette("matrix", {
        ids <- retrieveIds("ipi", "rna_seq")
        
        mat <- retrieveMatrix(
            "ipi",
            "rna_seq",
            recordNames = ids[grep("live",ids)[1:25]],
            attributeNames = "gene_counts")
    })
    
    expect_type(mat, "double")
    expect_equal(ncol(mat), 25)
    expect_gt(nrow(mat), 50000)
    expect_true(all(colnames(mat) %in% ids))
    expect_true(all(grepl("ENSG",rownames(mat)[1:100])))
})

