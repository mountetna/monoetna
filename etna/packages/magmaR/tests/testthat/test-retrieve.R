# This code tests the retrieve and retrieveJSON functions
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-retrieve.R")

test_that("retrieve & retrieveJSON work with minimal input", {
    vcr::use_cassette("retrieve_ALLs", {
        df <- retrieve(
            "ipi", "experiment")
    })
    
    expect_s3_class(df, "data.frame")
    expect_true(nrow(df) >= 20 && ncol(df) > 4)
    
    vcr::use_cassette("retrieveJSON_ALLs", {
        json <- retrieveJSON(
            "ipi", "experiment")
    })
    
    expect_type(json, "list")
    expect_equal(nrow(df), length(json$models$experiment$documents))
    expect_gt(length(json$models$experiment$documents[[1]]), 4)
})

test_that("retrieve works with targetted input", {
    vcr::use_cassette("retrieve_rec3_att1", {
        ret <- retrieve(
            "ipi", "experiment",
            recordNames = c("Adrenal", "Bladder", "Breast"),
            attributeNames = c("description"))
    })
    
    expect_s3_class(ret, "data.frame")
    # Id attribute is NOT targetted, +1
    expect_equal(dim(ret), c(3,2))
})

test_that("retrieve works with targetted input", {
    vcr::use_cassette("retrieve_rec1_att3", {
        ret <- retrieve(
            "ipi", "experiment",
            recordNames = c("Adrenal"),
            attributeNames = c("description", "name", "project"))
    })
    
    expect_s3_class(ret, "data.frame")
    # Id attribute is targetted
    expect_equal(dim(ret), c(1,3))
})
