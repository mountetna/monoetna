# This code tests the updateValues & updateMatrix functions
# library(magmaR); library(testthat); source("tests/testthat/helper-magmaR.R"); source("tests/testthat/test-retrieve.R")

### NOTE: updateValues() is currently very simple and is used by updateMatrix(), so no direct tests are needed for that fxn.

vcr::use_cassette("update_1", {
    test_that("updateMatrix can take in data directly", {
        
        mat <<- retrieveMatrix("example", "rna_seq", "all", "gene_counts")
        
        expect_output(
            x <- updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat,
                auto.proceed = TRUE),
            "This update() will update data for 12 records.\n/update: successful.", fixed = TRUE
        )
        mat_after <<- retrieveMatrix("example", "rna_seq", "all", "gene_counts")
    
        expect_identical(mat,mat_after)
    })
})

vcr::use_cassette("update_2", {
    test_that("updateMatrix can take in data as csv", {
        
        expect_output(
            x <- updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.csv",
                auto.proceed = TRUE),
            "This update() will update data for 12 records.\n/update: successful.", fixed = TRUE
        )
        mat_after <- retrieveMatrix("example", "rna_seq", "all", "gene_counts")
    
        expect_identical(mat,mat_after)
    })
})

vcr::use_cassette("update_3", {
    test_that("updateMatrix can take in data as tsv", {
        
        expect_output(
            x <- updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.tsv",
                separator = "\t",
                auto.proceed = TRUE),
            "This update() will update data for 12 records.\n/update: successful.", fixed = TRUE
        )
        mat_after <- retrieveMatrix("example", "rna_seq", "all", "gene_counts")
    
        expect_identical(mat,mat_after)
        
        # Error if 'separator' is not changed:
        expect_error(
            x <- updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.tsv",
                auto.proceed = TRUE),
            "Parsing error.", fixed = TRUE
        )
        
    })
})

vcr::use_cassette("update_4", {
    test_that("updateMatrix gives expected errors/warnings", {
        
        # Error when matrix does not have colnames (recordNames)
        mat_noIDs <- mat
        colnames(mat_noIDs) <- NULL
        expect_error(
            updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat_noIDs),
            "Colnames of matrix should be record names.", fixed = TRUE
        )
        
        # Error when matrix rownames (features / magma "colnames") do not match attribute 'options'
        mat_notGENEs <- mat
        rownames(mat_notGENEs) <- paste0("NOPE", seq_len(nrow(mat)))
        expect_error(
            updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat_notGENEs),
            "Validation error: rownames of 'matrix' are not valid 'options' for", fixed = TRUE
        )
        
        # Error when run non-interactively without 'auto.proceed = TRUE' or when user says "no" to proceeding.
        # Note: when running interactively, SAY NO!
        skip_if(interactive(), "running in interactive mode")
        expect_warning(
            updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat),
            "To run in non-interactive mode, set 'auto.proceed = TRUE'.", fixed = TRUE)
    })
})

### Note: suppressWarnings added below because auto.proceed was removed so that
# none of the calls will actually cause an /update to happen. The warning being
# suppressed is the one which tells users how to allow the function to work
# even in the non-interactive mode.

### Note2: When running in an interactive mode, always say NO to proceeding if asked.

vcr::use_cassette("update_5", {
    test_that("updateMatrix gives expected messages", {
        
        # When all records already exist:
        expect_output(
            suppressWarnings(
                updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat)),
"This update() will update data for 12 records.", fixed = TRUE
        )
        
        # When some records are new and some are not:
        mat_halfIDs_wrong <- mat
        colnames(mat_halfIDs_wrong)[1:6] <- paste0("WRONG", 1:6)
        expect_output(
            suppressWarnings(
                updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat_halfIDs_wrong)),
"This update() will create 6 new records: WRONG1, WRONG2, WRONG3, WRONG4, WRONG5, WRONG6 
WARNING: Check the above carefully. Once created, there is currently no way to remove records form magma.

This update() will update data for 6 records.", fixed = TRUE
        )
        
        # When user-cancels the upload (or in non-interactive mode).
        expect_equal(
            suppressWarnings(
                updateMatrix(projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat)),
            "No /update performed."
        )
        
        ###### Add test for message when attribute has no 'options'
        # skip_if(!identical(mat, mat_after), "before/after updateMatrix() are not equal")
    })
})

