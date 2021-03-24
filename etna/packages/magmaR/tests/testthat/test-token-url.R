# This code tests the magmaRset function
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-token-url.R")

test_that("token, url, and opts get stored in expected spots, & 'followlocation' added when left out by user", {
    
    expect_message(
        targ <- magmaRset(
            token = Sys.getenv("TOKEN"),
            url = Sys.getenv("URL"),
            opts = list(option1 = FALSE)),
        "'followlocation = FALSE' added",
        fixed = TRUE)
    
    expect_identical(
        targ$token,
        Sys.getenv("TOKEN"))
    expect_identical(
        targ$url,
        Sys.getenv("URL"))
    expect_identical(
        targ$opts,
        list(option1 = FALSE, followlocation = FALSE))
})

test_that("default, production, magma link is filled in when not gven", {
    
    expect_identical(
        magmaRset("")$url,
        "https://magma.ucsf.edu")
})
