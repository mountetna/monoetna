# This code tests functions of utils.R
# library(magmaR); library(testthat); source("R/utils.R"); source("tests/testthat/test-utils.R")

test_that("A simple revisions structure is converted to JSON in expected way", {
    # Single string
    expect_true(
        '{"model_name":{"record_name":{"attribute_name":"hello"}}}' == .jsonify(list(
            'model_name' = list(
                'record_name' = list(
                    'attribute_name' = 'hello'
                )
            )
        ))
    )
    # String vector
    expect_true(
        '{"model_name":{"record_name":{"attribute_name":["x","y"]}}}' == .jsonify(list(
            'model_name' = list(
                'record_name' = list(
                    'attribute_name' = c('x', 'y')
                )
            )
        ))
    )
    # Single number
    expect_true(
        '{"model_name":{"record_name":{"attribute_name":1}}}' == .jsonify(list(
            'model_name' = list(
                'record_name' = list(
                    'attribute_name' = 1
                )
            )
        ))
    )
    # Boolean
    expect_true(
        '{"model_name":{"record_name":{"attribute_name":true}}}' == .jsonify(list(
            'model_name' = list(
                'record_name' = list(
                    'attribute_name' = TRUE
                )
            )
        ))
    )
})

### Special Value Cases:

# NULL, special value in R which we want to map to the JSON null
test_that("NULL is converted to JSON in expected way, as null", {
    expect_true(
        '{"record_name":null}' == .jsonify(list(
            'record_name' = NULL
        ))
    )
})

# NA, a special value in R which we want to map to the string "NA" so that null
# remains reserved for blanking.
test_that("NA is converted to JSON in expected way, as 'NA' string", {
    expect_true(
        '{"record_name":"NA"}' == .jsonify(list(
            'record_name' = NA
        ))
    )
})

# [] for retrieveTemplate() special case 
test_that(" 'recordNames = []' of retrieve calls converted to JSON in expected way, as non-string []", {
    expect_true(
        '{"record_name":[]}' == .jsonify(list(
            'record_name' = .match_expected_recName_structure("[]")
        ))
    )
})

# "idenifier" versus actual attribute names for retrieve*() attributeNames
test_that(" 'attributeNames' special cases of retrieve converted to JSON in expected ways", {
    expect_true(
        '{"attribute_name":"identifier"}' == .jsonify(list(
            'attribute_name' = .match_expected_attName_structure("identifier")
        ))
    )
    expect_true(
        '{"attribute_name":["name"]}' == .jsonify(list(
            'attribute_name' = .match_expected_attName_structure("name")
        ))
    )
    expect_true(
        '{"attribute_name":["name","other"]}' == .jsonify(list(
            'attribute_name' = .match_expected_attName_structure(c("name", "other"))
        ))
    )
})
