# This code tests the retrieveMatrix function
# library(magmaR); library(testthat); source("tests/testthat/helper-magmaR.R"); source("tests/testthat/test-matrix.R")

targ <- magmaRset(
    token = TOKEN,
    url = URL)

vcr::use_cassette("matrix", {
    
    ids <- retrieveIds(targ, "example", "rna_seq")
    
    test_that("retrieveMatrix", {
    
        mat <- retrieveMatrix(
            targ,
            "example",
            "rna_seq",
            attributeNames = "gene_counts")
    
        expect_type(mat, "integer")
        expect_equal(dim(mat), c(40,12))
        
        # Column names = identifiers
        expect_true(all(colnames(mat) %in% ids))
        # Row names = gene names (pulled from template)
        expect_equivalent(rownames(mat), paste0("gene", 1:40))
    })

    test_that("retrieveMatrix warns but ignores empty records", {
        
        expect_warning(
            mat <- retrieveMatrix(
                targ,
                "example",
                "rna_seq",
                recordNames = c(ids, "not_a_record"),
                attributeNames = "gene_counts"),
            "Empty record, not_a_record, was ignored."
        )
        
        expect_type(mat, "integer")
        expect_equal(dim(mat), c(40,12))
    })
})
