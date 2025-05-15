# library(testthat); for (util_file in list.files("R", pattern="^util", full.names = TRUE)) { source(util_file) }; rm(util_file); source("tests/testthat/setup.R", chdir = TRUE); source("R/subsetting.R")

con_num1_null <- cons[[1]]
con_num2_and <- cons2[[4]]
con_str_or <- list(
    col = "cell_type",
    def = list("B", "Myeloid"),
    logic = "OR"
)
con_str_and <- con_str_or
con_str_and$logic <- "AND"
con_bool_or <- list(
    col = "has_thing",
    def = list("true"),
    logic = "OR"
)
con_bool_and <- list(
    col = "has_thing",
    def = list("true"),
    logic = "AND"
)

test_that(".check_constraint_format allows proper formats, and outputs expected types", {
    # numeric
    expect_equal(
        .check_constraint_format(con_num1_null, 1),
        "numeric"
    )
    expect_equal(
        .check_constraint_format(con_num2_and, 2),
        "numeric"
    )

    expect_equal(
        .check_constraint_format(con_str_or, 2),
        "string"
    )

    expect_equal(
        .check_constraint_format(con_bool_or, 2),
        "boolean"
    )
})

test_that(".check_constraint_format catches missing logic", {
    expect_error(
        .check_constraint_format(con_num1_null, 1),
        NA
    )
    expect_error(
        .check_constraint_format(con_num1_null, 2),
        "'logic' missing;"
    )
})

test_that(".interpret_contraint gives expected outputs", {
    expect_equal(
        which(.interpret_constraint(con_num1_null$col, con_num1_null$def, "numeric", df)),
        c(1,3)
    )
    expect_equal(
        which(.interpret_constraint(con_num2_and$col, con_num2_and$def, "numeric", df)),
        c(2:4,6)
    )
    expect_equal(
        which(.interpret_constraint(con_str_or$col, con_str_or$def, "string", df)),
        c(2,3,5,6)
    )
    expect_equal(
        which(.interpret_constraint(con_bool_or$col, con_bool_or$def, "boolean", df)),
        c(1:3)
    )
})

test_that("subset function gives expected outputs", {
    expect_equal(
        subsetDF_index_targets(df, list(con_num1_null)),
        c(1,3)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_num2_and)),
        c(2:4,6)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_str_or)),
        c(2,3,5,6)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_bool_or)),
        c(1:3)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_num1_null, con_str_or)),
        c(1:3,5,6)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_num1_null, con_str_and)),
        c(3)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_str_or, con_bool_and)),
        c(2,3)
    )
    expect_equal(
        subsetDF_index_targets(df, list(con_bool_and, con_str_or)),
        c(1:3,5,6)
    )
})
