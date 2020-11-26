# This code tests the query function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-query.R")

test_that("query_list", {
    vcr::use_cassette("query_list", {
        # query obtains sample ids for records with rna_seq data
        qry <- query(
            "ipi",
            queryTerms = 
                list('rna_seq',
                     list('::has', 'gene_tpm'),
                     '::all',
                     'sample',
                     '::identifier'
                ))
    })
    
    expect_type(qry, "list")
})

test_that("query_df", {
    vcr::use_cassette("query_df", {
        # query obtains sample ids for records with rna_seq data
        qry <- query(
            "ipi",
            queryTerms = 
                list('rna_seq',
                     list('::has', 'gene_tpm'),
                     '::all',
                     'sample',
                     '::identifier'
                ),
            format = "df")
    })
    
    vcr::use_cassette("sample_ids", {
        # query obtains sample ids for records with rna_seq data
        ids <- retrieveIds("ipi", "sample")
    })
    
    expect_s3_class(qry, "data.frame")
    expect_true(all(qry[,2] %in% c(ids, NA)))
})
