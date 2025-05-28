# library(testthat); for (util_file in list.files("R", pattern="^util", full.names = TRUE)) { source(util_file) }; rm(util_file); source("tests/testthat/setup.R", chdir = TRUE)

test_that(".is_checks work", {
    # .is_numeric
    expect_equal(
        .is_numeric(test_vector_na),
        c(rep(TRUE, 5),
          rep(FALSE, 3),
          rep(FALSE, 2),
          rep(FALSE, 4),
          TRUE)
    )
    expect_equal(
        .is_numeric(test_vector_na, FALSE),
        c(rep(TRUE, 5),
          rep(FALSE, 3),
          rep(FALSE, 2),
          rep(FALSE, 4),
          FALSE)
    )
    # .is_logical
    expect_equal(
        .is_logical(test_vector_na),
        c(rep(FALSE, 5),
          rep(FALSE, 3),
          FALSE, TRUE,
          rep(TRUE, 4),
          TRUE)
    )
    expect_equal(
        .is_logical(test_vector_na, FALSE),
        c(rep(FALSE, 5),
          rep(FALSE, 3),
          FALSE, TRUE,
          rep(TRUE, 4),
          FALSE)
    )
})
