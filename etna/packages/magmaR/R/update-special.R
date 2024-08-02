#' A matrix-specific wrapper of \code{\link{updateValues}}
#' @description A matrix-specific wrapper of \code{\link{updateValues}} which can take in a matrix, data.frame, or file path, directly.
#' @inheritParams retrieve
#' @inheritParams updateValues
#' @param attributeName String naming the matrix attribute for which to upload data.
#' @param matrix A matrix or dataframe containing the data to upload to magma.
#' 
#' Alternatively, a String specifying the file path of a file containing such data.
#' 
#' No matter the provision method, colnames must be record identifiers, and rownames should match the values of 'options' associated with the target 'attribute'.
#' Check the 'See Also' section below for how to determine the needed 'options'.
#' @param separator String indicating the field separator to use if providing \code{matrix} as a file path.
#' Default = \code{","}.
#' @param revisions.only Logical. For troubleshooting purposes, when set to \code{TRUE}, no data will be sent to magma.
#' Instead, the list structure that would have been passed to the \code{revisions} input of \code{\link{updateValues}} is returned as output.
#' @return None directly.
#' 
#' The function sends data to magma, and the only outputs are information reported via the console.
#' @details This function utilizes the magma/query function, documented here \url{https://mountetna.github.io/magma.html#update},
#' to upload data to a matrix attribute (named \code{attributeName}) of the \code{modelName} model of \code{projectName} project.
#' 
#' \code{matrix} data are provided either as a matrix, dataframe, or file path which points toward such data.
#' If given as a file path, the \code{separator} input can be used to adjust for whether the file is a csv (the default, \code{separator = ","}), or tsv, \code{separator = "\t"}, or other format.
#' 
#' Data is then validated by ensuring that all row names are among the valid 'options' of the target attribute (See the See Also section below for a note on how to explore these options yourself.).
#' Rows are reordered to be in the same order as these 'options'.
#' 
#' For any missing 'options', rows of NAs are added.
#' 
#' The data is then transformed and passed along to \code{\link{updateValues}}.
#' 
#' The \code{updateValues()} function will summarize records to be updated and allow the user to double-check this information before proceeding.
#' 
#' This user-prompt step can be bypassed (useful when running in a non-interactive way) by setting \code{auto.proceed = TRUE}, but NOTE:
#' It is a good idea to always check carefully before proceeding, if possible.
#' Data can be overwritten with NAs or zeros or the like, but improperly named records cannot be easily removed.
#' 
#' @seealso
#' \code{\link{updateFromDF}} for a more flexible function for uploading multiple attributes-worth of (non-matrix) data at a time.
#' 
#' \code{\link{updateValues}} for the more direct replica of \code{magma/update} which is more even more flexible that \code{updateFromDF}, though a bit more complicated to use.
#' 
#' \code{\link{retrieveTemplate}}, then check the \code{<output>$models$<modelName>$template$attributes$<attributeName>$options} to explore the rownames that your matrix should have.
#'
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying \code{magma/update} function.
#' 
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     ### Note that you likely do not have write-permissions for the 'example'
#'     # project, so this code can be expected to give an authorization error.
#'     
#'     ### Retrieve some data from magma, then update that same data.
#'     mat <- retrieveMatrix(magma, "example", "rna_seq", "all", "gene_tpm")
#' 
#'     updateMatrix(
#'         target = magma,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         attributeName = "gene_tpm",
#'         matrix = mat)
#' }
#'
#' @importFrom utils read.csv
updateMatrix <- function(
    target,
    projectName,
    modelName,
    attributeName,
    matrix,
    autolink = FALSE,
    dryRun = FALSE,
    separator = ",",
    auto.proceed = FALSE,
    revisions.only = FALSE,
    template = NULL,
    ...) {
    
    ### Read in matrix if passed as a string (file-location)
    if (is.character(matrix)) {
        matrix <- read.csv(
            matrix, sep = separator, row.names = 1, check.names = FALSE)
        if (ncol(matrix)<1){
            stop("Parsing error. Is 'separator' correct?")
        }
    }
    # Validate that there are column names
    if (is.null(colnames(matrix))) {
        stop("Colnames of matrix should be record names.")
    }
    
    ### Validate rownames (genes / value names)
    # Obtain 'validation' / rowname options
    if (identical(template, NULL)) {
        template <- retrieveTemplate(target, projectName)
    }
    row_options <- 
        template$models[[modelName]]$template$attributes[[attributeName]]$validation$value
    
    if (length(row_options)<1) {
        cat("WARNING: Target attribute does not have 'validation' info: no feature-names validation can be performed.\n\n")
    } else {
        # Check rownames against options
        if (! all(rownames(matrix) %in% row_options)) {
            stop("Validation error: rownames of 'matrix' are not valid options for ", attributeName,".")
        } else {
            # Add NAs for any options not in the current matrix
            for (not_there in setdiff(row_options, rownames(matrix))) {
                matrix <- rbind(matrix, NA)
                rownames(matrix)[nrow(matrix)] <- not_there
            }
            # Reorder as in options
            matrix <- matrix[row_options, ,drop = FALSE]
        }
    }
    
    # Transform into a nested list
    # Note: Do not supply recordNames directly to vapply as any "-" will be
    #   converted to "."
    rec_att_vals <- vapply(seq_len(ncol(matrix)), function(x) {
        x <- list(matrix[,x])
        names(x) <- attributeName
        list(x)
        }, FUN.VALUE = list(1))
    names(rec_att_vals) <- colnames(matrix)
    
    revs <- list(rec_att_vals)
    # Because list(modelName = ...) would not substitute the value of modelName
    names(revs) <- modelName
    
    if (revisions.only) {
        return(revs)
    }
    ### Pass to updateValues() to:
    # Summarize to user
    # Check with the user before proceeding
    # Perform upload
    updateValues(
        target = target,
        projectName = projectName,
        revisions = revs,
        auto.proceed = auto.proceed,
        autolink = autolink,
        dryRun = dryRun,
        template = template,
        ...)
}

#' Easier to use wrapper of \code{\link{updateValues}}
#' @description A wrapper of \code{\link{updateValues}} which takes in updates in the form of a dataframe, csv, tsv, with rows = records and columns = attributes.
#' @inheritParams retrieve
#' @inheritParams updateMatrix
#' @inheritParams updateValues
#' @param df A dataframe, containing the data to upload to magma.
#' 
#' Alternatively, a String specifying the file path of a file containing such data.
#' 
#' See below for additional formatting details.
#' @param separator String indicating the field separator to use if providing \code{df} as a file path.
#' Default = \code{","}.
#' Use \code{"\t"} for tsvs.
#' @param table.method "replace" or "append".
#' Sets the methodology used for building the \code{revisions} to request for table-type model updates:
#' \itemize{
#' \item "append": Add the currently defined records, attaching them to parents defined in the \code{df}.
#' \item "replace": Add the currently defined records, attaching them to parents defined in the \code{df}, AND unlink / remove all current records, of \code{modelName}, attached to parent records defined in the \code{df}.
#' }
#' @return None directly.
#' 
#' The function sends data to magma, and the only outputs are information reported via the console.
#' @details This function provides a simple method for updating multiple attributes of multiple magma records provided as a rectangular dataframe, or equivalent file structure.
#' It utilizes the magma/query function, documented here \url{https://mountetna.github.io/magma.html#update},
#' to upload data after converting to the format required by that function.
#' 
#' The user-indicated \code{df} is read in, presented to the user for inspection, then transformed to the necessary format and passed along to \code{\link{updateValues}}.
#' 
#' The \code{updateValues()} function will then summarize records to be updated and allow the user to double-check this information before proceeding.
#' 
#' This user-prompt step can be bypassed (useful when running in a non-interactive way) by setting \code{auto.proceed = TRUE}, but NOTE:
#' It is a good idea to always check carefully before proceeding, if possible.
#' Data can be overwritten with NAs or zeros or the like, or disconnected from parent records, but improperly named records cannot be easily removed.
#' 
#' For "standard" models with explicit identifiers, the function targets the \code{df}'s row-indicated records and column-indicated attributes of the \code{modelName} model of \code{projectName} project.
#' In such cases, the first column of \code{df} must contain the identifiers of the records your wish to update.
#' 
#' For table-type models which do not have explicit identifiers, the function creates records per each row of the given \code{df} with the requested values filled in for column-indicated attributes of the \code{modelName} model of \code{projectName} project, and attaches these records to the indicated parent records.
#' In such cases, a column named as the parent model must exist in \code{df} to provide the parent identifiers of all requested new data.
#' In such cases, \code{table.method} must also be given as either \code{"append"} or \code{"replace"}.
#' 
#' \code{df} can be provided either as a data.frame directly, or as a file path pointing to a file containing such data.
#' If given as a file path, the \code{separator} input can be used to adjust for whether the file is a csv (the default, \code{separator = ","}), or tsv, \code{separator = "\t"}, or other format.
#' 
#' The \code{df} data structure when targeting 'standard' models:
#' \itemize{
#' \item Rows = records, with the first column indicating the record identifiers.
#' \item Columns = represent the data desired to be given for each attribute.
#' \item Column Names (or the top row when providing a file) = attribute names.
#' Except for the first column (ignored as this column's data are used as identifiers), all column names must be valid attribute names of the target \code{modelName}.
#' }
#' 
#' The \code{df} data structure when targeting table-type models:
#' \itemize{
#' \item Rows = records, but no identifiers are needed.
#' \item Columns = represent the data desired to be given for each attribute.
#' \item Column Names (or the top row when providing a file) = attribute names.
#' At least one column must be named after the parent model and must represent parent model identifiers.
#' }
#' 
#' @section Use Case. Using this function to change records' identifiers:
#' 
#' To do so, provide a file or dataframe where
#' 1) The first column, named something random Iits name will be ignored.), contains current identifiers;
#' 2) Some other column, named as the attribute which is treated as the identifier for the model, contains the new identifiers
#' 
#' To determine the identifier attribute's name, you can use \code{\link{retrieveTemplate}}:
#' 
#' \code{retrieveTemplate(<target>, <projectName>)$models$<modelName>$template$identifier}.
#' 
#' 
#' @seealso
#' \code{\link{updateMatrix}} for uploading matrix data
#' 
#' \code{\link{updateValues}} for a more direct replica of \code{magma/update} which is more flexible, though a bit more complicated to use.
#' 
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying \code{magma/update} function.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     ### Note that you likely do not have write-permissions for the 'example'
#'     # project, so this code can be expected to give an authorization error.
#'     
#'     ### Retrieve some data from magma, which will be in the proper format.
#'     df <- retrieve(
#'         magma, projectName = "example", modelName = "rna_seq",
#'         recordNames = "all",
#'         attributeNames = c("tube_name", "biospecimen", "cell_number")
#'         )
#'     df
#'     
#'     updateFromDF(
#'         target = magma,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         df = df)
#' }
#'
#' @importFrom utils read.csv
#' @importFrom stats setNames
updateFromDF <- function(
    target,
    projectName,
    modelName,
    df,
    table.method = NULL,
    autolink = FALSE,
    dryRun = FALSE,
    separator = ",",
    auto.proceed = FALSE,
    revisions.only = FALSE,
    template = NULL,
    ...) {
    
    ### Read in df if passed as a string (file-location)
    if (is.character(df)) {
        df <- read.csv(
            df, sep = separator, check.names = FALSE, header = TRUE)
        if (ncol(df)<2){
            stop("Parsing error. Is 'separator' correct?")
        }
    }
    cat("Data recieved:\n")
    print(df)
    
    if (identical(template, NULL)) {
        template <- retrieveTemplate(target, projectName)
    }
    isTable <- .is_table_model(target, projectName, modelName, template)
    
    ### Validation & convert revisions
    if (isTable) {
        ### Targeting table models, no identifiers and special replace case
        if (identical(table.method, NULL)) {
            stop("'", modelName, "' model of '", projectName, "' project is a table-type model, but 'table.method' input is not set. Set it to \"append\" or \"replace\".")
        }
        if (!table.method %in% c("append", "replace")) {
            stop("'table.method' must be \"append\" or \"replace\", but is: \"", table.method, "\".")
        }
        parentModelName <- template$models[[modelName]]$template$parent
        if (!parentModelName %in% names(df)) {
            stop("Parent attribute, ", parentModelName, ", must be given when updating records of table-type models.")
        }

        # Prep revisions
        tempIDs <- paste0("::temp-id-", seq_len(nrow(df)))
        df <- cbind('__IDS__' = tempIDs, df)
        revs <- .df_to_revisions(df, modelName = modelName)

        # Add to revisions for 'replace' method
        if (table.method=="replace") {
            parents <- unique(df[[parentModelName]])
            parent_revs <- lapply(
                parents,
                function(this) {
                    setNames(list(tempIDs[df[[parentModelName]]==this]), modelName)
                }
            )
            names(parent_revs) <- parents
            revs[[parentModelName]] <- parent_revs
        }
    } else {
        ### Targeting "standard" models, with identifiers
        # Validate that 1st column, which should be IDs, is all unique values
        if (any(duplicated(df[,1]))) {
            stop("Values of 1st column (record identifiers) must be unique.")
        }
        if (ncol(df) < 2) {
            stop("df has one column, but at least 2 are required.")
        }

        # Prep revisions
        revs <- .df_to_revisions(df, modelName = modelName)
    }
    
    if (revisions.only) {
        return(revs)
    }
    ### Pass to updateValues() to:
    # Summarize to user
    # Check with the user before proceeding
    # Perform upload
    updateValues(
        target = target,
        projectName = projectName,
        revisions = revs,
        autolink = autolink,
        dryRun = dryRun,
        auto.proceed = auto.proceed,
        template = template,
        ...)
}
