# This code tests the retrieveModels, retrieveIds, retrieveAttributes, and retrieveTemplate functions
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-special.R")

test_that("retrieveTemplate", {
  vcr::use_cassette("Temp", {
      ret <- retrieveTemplate("example")
  })
  expect_type(ret, "list")
})

test_that("retrieveModels", {
  vcr::use_cassette("models", {
      ret <- retrieveModels("example")
  })
  expect_type(ret, "character")
})

test_that("retrieveIds", {
  vcr::use_cassette("ids", {
      ret <- retrieveIds("example", "subject")
  })
  expect_type(ret, "character")
})

test_that("retrieveAttributes", {
  vcr::use_cassette("Atts", {
      ret <- retrieveAttributes("example", "subject")
  })
  expect_type(ret, "character")
})
