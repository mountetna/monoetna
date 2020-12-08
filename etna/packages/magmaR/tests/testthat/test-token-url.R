# This code tests the .get_TOKEN & .get_URL functions
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-token-url.R")

test_that("automated token setting gives error in non-interactive mode", {
    
    rm(.MAGMAR_TOKEN, envir = .GlobalEnv)
    skip_if(interactive())
    
    expect_error(
        magmaR:::.get_TOKEN(),
        "Please provide", fixed = TRUE)
})

test_that("automated token setting retrieves .MAGMAR_TOKEN", {
    
    .GlobalEnv$.MAGMAR_TOKEN <- Sys.getenv("TOKEN")
    
    expect_identical(
        magmaR:::.get_TOKEN(),
        Sys.getenv("TOKEN"))
})

test_that("automated token setting creates & retrieves .MAGMAR_TOKEN", {
    rm(.MAGMAR_URL, envir = .GlobalEnv)
    
    expect_identical(
        magmaR:::.get_URL(),
        "https://magma.ucsf.edu")
    expect_identical(
        .MAGMAR_URL,
        "https://magma.ucsf.edu")
})
