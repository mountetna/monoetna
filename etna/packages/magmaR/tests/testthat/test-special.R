# This code tests the retrieveModels, retrieveIds, retrieveAttributes, and retrieveTemplate functions
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-special.R")

test_that("retrieveTemplate", {
  vcr::use_cassette("Temp", {
      ret <- retrieveTemplate("ipi")
  })
  expect_type(ret, "list")
})

test_that("retrieveModels", {
  vcr::use_cassette("models", {
      ret <- retrieveModels("ipi")
  })
  expect_type(ret, "character")
})

test_that("retrieveIds", {
  vcr::use_cassette("ids", {
      ret <- retrieveIds("ipi", "experiment")
  })
  expect_type(ret, "character")
})

test_that("retrieveAttributes", {
  vcr::use_cassette("Atts", {
      ret <- retrieveAttributes("ipi", "experiment")
  })
  expect_type(ret, "character")
})
