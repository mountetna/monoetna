# This code tests the retrieveModels, retrieveIds, retrieveAttributes, and retrieveTemplate functions
# library(magmaR); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-special.R")

targ <- magmaRset(
    token = TOKEN,
    url = URL)

test_that("retrieveTemplate", {
  vcr::use_cassette("Temp", {
      ret <- retrieveTemplate(targ, "example")
  })
  expect_type(ret, "list")
})

test_that("retrieveModels", {
  vcr::use_cassette("models", {
      ret <- retrieveModels(targ, "example")
  })
  expect_type(ret, "character")
})

test_that("retrieveIds", {
  vcr::use_cassette("ids", {
      ret <- retrieveIds(targ, "example", "subject")
  })
  expect_type(ret, "character")
})

test_that("retrieveAttributes", {
  vcr::use_cassette("Atts", {
      ret <- retrieveAttributes(targ, "example", "subject")
  })
  expect_type(ret, "character")
})
