# This code tests the retrieveMatrix function
# library(magmaR); library(testthat); source("tests/testthat/helper-magmaR.R"); source("tests/testthat/test-matrix.R")

vcr::use_cassette("ids_rna", {
    ids <- retrieveIds("ipi", "rna_seq")
})

test_that("retrieveMatrix", {
    vcr::use_cassette("matrix", {
        mat <- retrieveMatrix(
            "ipi",
            "rna_seq",
            recordNames = ids[grep("live",ids)[1:12]],
            attributeNames = "gene_counts")
    })
    
    expect_type(mat, "double")
    expect_equal(ncol(mat), 12)
    expect_gt(nrow(mat), 50000)
    expect_true(all(colnames(mat) %in% ids))
    expect_true(all(grepl("ENSG",rownames(mat)[1:100])))
})

vcr::use_cassette("matrix_2", {
    test_that("retrieveMatrix warns but ignores empty records", {
        expect_warning(
            mat <- retrieveMatrix(
                "ipi",
                "rna_seq",
                recordNames = c(ids[grep("live",ids)[1:2]], "not_a_record"),
                attributeNames = "gene_counts"),
            "Empty record, not_a_record, was ignored."
        )
        
        expect_type(mat, "double")
        expect_equal(ncol(mat), 2)
        expect_gt(nrow(mat), 50000)
        expect_true(all(colnames(mat) %in% ids))
        expect_true(all(grepl("ENSG",rownames(mat)[1:100])))
    })
})
