# This code tests the query function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-query.R")

test_that("query_list", {
    vcr::use_cassette("query_list", {
        # query obtains sample ids for records with rna_seq data
        qry <- query(
            "example",
            queryTerms = 
                list('rna_seq',
                     list('::has', 'gene_tpm'),
                     '::all',
                     'subject',
                     '::identifier'
                ))
    })
    
    expect_type(qry, "list")
})

test_that("query_df", {
    vcr::use_cassette("query_df", {
        
        ids <- retrieveIds("example", "subject")
        
        # query obtains sample ids for records with rna_seq data
        qry <- query(
            "example",
            queryTerms = 
                list('rna_seq',
                     list('::has', 'gene_tpm'),
                     '::all',
                     'subject',
                     '::identifier'
                ),
            format = "df")
    })
    
    expect_s3_class(qry, "data.frame")
    expect_true(all(qry[,2] %in% c(ids, NA)))
})
