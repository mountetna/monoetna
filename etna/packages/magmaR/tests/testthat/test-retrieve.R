# This code tests the retrieve and retrieveJSON functions
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-retrieve.R")

.GlobalEnv$.MAGMAR_TOKEN <- Sys.getenv("TOKEN")

test_that("retrieve works with minimal input", {
    vcr::use_cassette("retrieve_targetted", {
        ret <- retrieve(
            "ipi", "experiment",
            recordNames = c("Adrenal", "Bladder", "Breast"),
            attributeNames = c("description"))
    })
    
    expect_s3_class(ret, "data.frame")
    expect_equal(dim(ret), c(3,2))
})
