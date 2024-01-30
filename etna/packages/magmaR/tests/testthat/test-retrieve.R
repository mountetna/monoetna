# This code tests the retrieve and retrieveJSON functions
# library(magmaR); library(testthat); source("tests/testthat/helper-magmaR.R"); source("tests/testthat/test-retrieve.R")

targ <- magmaRset(
    token = TOKEN,
    url = URL)

test_that("retrieve & retrieveJSON work with minimal input", {
    vcr::use_cassette("retrieve_ALLs", {
        json <- retrieveJSON(
            targ, "example", "subject")
        
        df <- retrieve(
            targ, "example", "subject")
    })
    
    expect_s3_class(df, "data.frame")
    expect_type(json, "list")
    
    # No. records same for both
    expect_equal(nrow(df), length(json$models$subject$documents))
    
    # exact No. of records & No. Attributes is many per record for the df method
    expect_true(nrow(df) == 12 && ncol(df) >= 3)
    
    # Attributes per record is many per record for the json method
    expect_gt(length(json$models$subject$documents[[1]]), 3)
    
    # Template returned for json method
    expect_true("template" %in% names(json[[1]][[1]]))
})

test_that("retrieve works with targetted input, 1att", {
    vcr::use_cassette("retrieve_rec3_att1", {
        ret <- retrieve(
            targ, "example", "subject",
            recordNames = c("EXAMPLE-HS1", "EXAMPLE-HS2"),
            attributeNames = c("group"))
    })
    
    expect_s3_class(ret, "data.frame")
    
    # Proper numbers retrieved?
    # Id attribute is NOT targeted, +1 column
    expect_equal(dim(ret), c(2,2))
})

test_that("retrieve works with targetted input, 1rec", {
    vcr::use_cassette("retrieve_rec1_att3", {
        ret <- retrieve(
            targ, "example", "subject",
            recordNames = c("EXAMPLE-HS1"),
            attributeNames = c("biospecimen", "name", "group"))
    })
    
    expect_s3_class(ret, "data.frame")
    
    # Proper numbers retrieved?
    # Id attribute is targeted, no +1 to columns
    expect_equal(dim(ret), c(1,3))
})

test_that("retrieveJSON can hide templates", {
    vcr::use_cassette("retrieve_json_no_temp", {
        ret <- retrieveJSON(
            targ, "example", "subject",
            recordNames = c("EXAMPLE-HS1"),
            attributeNames = c("biospecimen"),
            hideTemplate = TRUE)
    })
    
    expect_type(ret, "list")
    expect_false("template" %in% names(ret[[1]][[1]]))
})

test_that("retrieve can pass show_disconnected", {
    default <- retrieveJSON(
        targ, "example", "subject",
        recordNames = c("EXAMPLE-HS1"),
        attributeNames = c("biospecimen"),
        json.params.only = TRUE)
    json <- retrieveJSON(
        targ, "example", "subject",
        recordNames = c("EXAMPLE-HS1"),
        attributeNames = c("biospecimen"),
        showDisconnected = TRUE,
        json.params.only = TRUE)
    df <- retrieve(
        targ, "example", "subject",
        recordNames = c("EXAMPLE-HS1"),
        attributeNames = c("biospecimen"),
        showDisconnected = TRUE,
        json.params.only = TRUE)

    expect_equal(
        !default$show_disconnected,
        json$show_disconnected)
    expect_equal(
        json$show_disconnected,
        df$show_disconnected)
})
