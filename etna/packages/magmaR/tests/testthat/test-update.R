# This code tests the updateValues & updateMatrix functions
# library(magmaR); library(testthat); source("tests/testthat/helper-magmaR.R"); source("tests/testthat/test-retrieve.R")

### NOTE: updateValues() is currently used by updateMatrix(), so no direct tests are needed for that fxn.

targ <- magmaRset(
    token = TOKEN,
    url = URL)

vcr::use_cassette("update_1", {
    mat <- retrieveMatrix(targ, "example", "rna_seq", "all", "gene_counts")
    
    test_that("updateMatrix can take in data directly, as csv, or as tsv", {
        
        # Run for direct, and check that nothing actually changed 
        expect_output(
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat,
                auto.proceed = TRUE),
            "/update: successful.",
            fixed = TRUE
        )
        mat_after <- retrieveMatrix(targ, "example", "rna_seq", "all", "gene_counts")
        expect_identical(mat,mat_after)
        
        # csv, tsv identical
        revs_mat_csv <-
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.csv",
                revisions.only = TRUE)
        revs_mat_tsv <-
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.tsv",
                separator = "\t",
                revisions.only = TRUE)
        expect_identical(revs_mat_csv, revs_mat_tsv)
        
        # Errors if 'separator' is not changed for tsv
        expect_error(
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.tsv"),
            "Parsing error.",
            fixed = TRUE
        )
        
        # Run one, and check that nothing actually changed 
        expect_output(
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = "rna_seq_counts.csv",
                auto.proceed = TRUE),
            "/update: successful.",
            fixed = TRUE
        )
    })
})

vcr::use_cassette("update_2", {
    df <- retrieve(targ, "example", "rna_seq", "all",
                   c("tube_name", "cell_number", "fraction"))
    
    test_that("updateFromDF can take in data directly, as csv, or as tsv", {
        revs_df_direct <-
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df,
                revisions.only = TRUE)
        revs_df_csv <-
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = "rna_seq_df.csv",
                revisions.only = TRUE)
        revs_df_tsv <-
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = "rna_seq_df.tsv",
                separator = "\t",
                revisions.only = TRUE)
        expect_identical(revs_df_direct, revs_df_csv)
        expect_identical(revs_df_csv, revs_df_tsv)
        
        # Column#1 / identifiers column name is ignored
        revs_df_csv_idName <-
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = "rna_seq_df_idName.csv",
                revisions.only = TRUE)
        expect_identical(revs_df_csv, revs_df_csv_idName)
        
        # Errors if 'separator' is not changed for tsv
        expect_error(
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                df = "rna_seq_df.tsv"),
            "Parsing error.",
            fixed = TRUE
        )
        
        # Run one, and check that nothing actually changed 
        expect_output(
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df,
                auto.proceed = TRUE),
            "For model \"rna_seq\", this update() will update 12 records",
            fixed = TRUE
        )
        df_after <- retrieve(targ, "example", "rna_seq", "all",
                             c("tube_name", "cell_number", "fraction"))
        expect_identical(df,df_after)
    })
})

vcr::use_cassette("update_3", {
    test_that("updateMatrix gives expected errors/warnings", {
        
        # Error when matrix does not have colnames (recordNames)
        mat_noIDs <- mat
        colnames(mat_noIDs) <- NULL
        expect_error(
            updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                         matrix = mat_noIDs),
            "Colnames of matrix should be record names.", fixed = TRUE
        )
        
        # Error when matrix rownames (features / magma "colnames") do not match attribute 'validation'
        mat_notGENEs <- mat
        rownames(mat_notGENEs) <- paste0("NOPE", seq_len(nrow(mat)))
        expect_error(
            updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                         matrix = mat_notGENEs),
            "Validation error: rownames of 'matrix' are not valid options for", fixed = TRUE
        )
        
        # Error when run non-interactively without 'auto.proceed = TRUE'.
        # Note: when running interactively, SAY NO!
        skip_if(interactive(), "running in interactive mode")
        expect_warning(
            updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                         matrix = mat),
            "To run in non-interactive mode, set 'auto.proceed = TRUE'.", fixed = TRUE)
    })
})

vcr::use_cassette("update_4", {
    test_that("updateFromDF gives expected errors/warnings", {
        
        # Error when matrix does not have colnames (recordNames).
        df_dupID <- df
        df_dupID[2,1] <- df[1,1]
        expect_error(
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df_dupID),
            "Values of 1st column (record identifiers) must be unique.", fixed = TRUE
        )
        
        # Error when no attribute columns are given.
        df_IDs_only <- df[,1,drop=FALSE]
        expect_error(
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df_IDs_only),
            "df has one column, but at least 2 are required.", fixed = TRUE
        )
        
        # Error when run non-interactively without 'auto.proceed = TRUE'.
        skip_if(interactive(), "running in interactive mode")
        expect_warning(
            updateFromDF(
                target = targ, projectName = "example",
                modelName = "rna_seq",
                df = df),
            "To run in non-interactive mode, set 'auto.proceed = TRUE'.", fixed = TRUE)
    })
})

### Note: suppressWarnings added below because auto.proceed was removed so that
# none of the calls will actually cause an /update to happen. The warning being
# suppressed is the one which tells users how to allow the function to work
# even in the non-interactive mode.

### Note2: When running in an interactive mode, always say NO to proceeding if asked.

vcr::use_cassette("update_5", {
    test_that("updates give expected summary / not-run messages", {
        
        # When all records already exist:
        expect_output(
            suppressWarnings(
                updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat)),
"For model \"rna_seq\", this update() will update 12 records:
    EXAMPLE-HS10-WB1-RSQ1
    EXAMPLE-HS11-WB1-RSQ1
    EXAMPLE-HS12-WB1-RSQ1
    EXAMPLE-HS1-WB1-RSQ1
    EXAMPLE-HS2-WB1-RSQ1
    EXAMPLE-HS3-WB1-RSQ1
    EXAMPLE-HS4-WB1-RSQ1
    EXAMPLE-HS5-WB1-RSQ1
    EXAMPLE-HS6-WB1-RSQ1
    EXAMPLE-HS7-WB1-RSQ1
    EXAMPLE-HS8-WB1-RSQ1
    EXAMPLE-HS9-WB1-RSQ1",
            fixed = TRUE
        )
        
        # When some records are new and some are not:
        mat_halfIDs_wrong <- mat
        colnames(mat_halfIDs_wrong)[1:6] <- paste0("WRONG", 1:6)
        expect_output(
            suppressWarnings(
                updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat_halfIDs_wrong)),
"For model \"rna_seq\", this update() will create (or update) 6 NEW (or orphan) records:
    WRONG1
    WRONG2
    WRONG3
    WRONG4
    WRONG5
    WRONG6
WARNING: Check the above carefully. Once created, there is no easy way to remove records from magma.
For model \"rna_seq\", this update() will update 6 records:
    EXAMPLE-HS4-WB1-RSQ1
    EXAMPLE-HS5-WB1-RSQ1
    EXAMPLE-HS6-WB1-RSQ1
    EXAMPLE-HS7-WB1-RSQ1
    EXAMPLE-HS8-WB1-RSQ1
    EXAMPLE-HS9-WB1-RSQ1", fixed = TRUE
        )
        
        # When 0 records are current:
        mat_allIDs_wrong <- mat
        colnames(mat_allIDs_wrong) <- paste0("WRONG", seq_len(ncol(mat)))
        expect_output(
            suppressWarnings(
                updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat_allIDs_wrong)),
"For model \"rna_seq\", this update() will create (or update) 12 NEW (or orphan) records:
    WRONG1
    WRONG2
    WRONG3
    WRONG4
    WRONG5
    WRONG6
    WRONG7
    WRONG8
    WRONG9
    WRONG10
    WRONG11
    WRONG12
WARNING: Check the above carefully. Once created, there is no easy way to remove records from magma.
For model \"rna_seq\", this update() will update 0 records.", fixed = TRUE
        )
        
        # When user-cancels the upload (or in non-interactive mode).
        expect_equal(
            suppressWarnings(
                updateMatrix(target = targ, projectName = "example", modelName = "rna_seq", attributeName = "gene_counts",
                    matrix = mat)),
            "No /update performed."
        )
        
        # When the attribute does not have 'options'
        expect_output(
            suppressWarnings(
                updateMatrix(target = targ, projectName = "example", modelName = "rna_seq",
                    attributeName = "fraction",
                    matrix = mat)),
"WARNING: Target attribute does not have 'validation' info: no feature-names validation can be performed.\n\n", fixed = TRUE
        )
    })
})

vcr::use_cassette("update_6", {
    mat <- retrieveMatrix(targ, "example", "rna_seq", "all", "gene_counts")
    df <- retrieve(targ, "example", "rna_seq", "all",
                   c("tube_name", "cell_number", "fraction"))

    test_that("update functions pass through autolink and dryRun params", {

        expect_equal(
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat, auto.proceed = TRUE,
                json.params.only = TRUE,
                autolink=FALSE)$autolink,
            !updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat, auto.proceed = TRUE,
                json.params.only = TRUE,
                autolink=TRUE)$autolink
        )
        expect_equal(
            updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat, auto.proceed = TRUE,
                json.params.only = TRUE,
                dryRun=FALSE)$dry_run,
            !updateMatrix(
                targ, projectName = "example",
                modelName = "rna_seq", attributeName = "gene_counts",
                matrix = mat, auto.proceed = TRUE,
                json.params.only = TRUE,
                dryRun=TRUE)$dry_run
        )

        expect_equal(
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df, auto.proceed = TRUE,
                json.params.only = TRUE,
                autolink=FALSE)$autolink,
            !updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df, auto.proceed = TRUE,
                json.params.only = TRUE,
                autolink=TRUE)$autolink
        )
        expect_equal(
            updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df, auto.proceed = TRUE,
                json.params.only = TRUE,
                dryRun=FALSE)$dry_run,
            !updateFromDF(
                targ, projectName = "example",
                modelName = "rna_seq",
                df = df, auto.proceed = TRUE,
                json.params.only = TRUE,
                dryRun=TRUE)$dry_run
        )
    })
})
