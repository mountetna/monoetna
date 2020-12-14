#' Analogous to the '/update' function of magma
#' @inheritParams updateMatrix
#' @param revisions A list of named lists containing the data to be updated.
#' 
#' List structure:
#' \itemize{
#' \item top level name(s): modelNames, can be 1 or more.
#' \item 2nd level name(s): recordNames, can be 1 or more.
#' \item 3rd level name(s) & contents: the attributes to update & the values to use.
#' }
#' 
#' See \url{https://mountetna.github.io/magma.html#update} for additional formatting details.
#' 
#' @return None directly.  The function sends data to magma, and the only output is whether that send was successful when \code{verbose} is set to \code{TRUE}.
#' @details This function mimics the activity of the magma/update function, documented here \url{https://mountetna.github.io/magma.html#update},
#' with the main difference being that the \code{revisions} input should be in nested-list format.
#' Internally, this function does little more than directly pass its inputs along to magma/update via a curl request.
#' 
#' @seealso
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying magma/update function.
#' 
#' \code{\link{updateMatrix}} for a matrix-dedicated version of this function which can be provided a matrix, or matrix's file location, directly.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     
#'     # Note that you likely do not have write-permissions for the 'example'
#'     #   project, so this code can be expected to give and authorization error.
#' 
#'     updateValues(
#'         projectName = "example",
#'         revisions = list(
#'             # model
#'             'rna_seq' = list(
#'                 # record
#'                 'EXAMPLE-HS1-WB1-RSQ1' = list(
#'                     # attribute
#'                     'fraction' = list(
#'                         # value(s)
#'                         "Tcells"
#'                         )
#'                     )
#'                 )
#'             )
#'         )
#' }
#'
updateValues <- function(
    projectName,
    revisions = list(),
    token = .get_TOKEN(),
    auto.proceed = FALSE,
    ...
) {
    
    ### Summarize per model
    lapply(names(revisions), function(model) {
        
        current_ids <- retrieveIds(
            projectName, model, token, ...)
        
        .summarize_model_values(revisions[[model]], model, current_ids)
    })
    
    ### Check if should move forward
    go <- .ask_before_proceeding(auto.proceed)
    if (! go %in% c("Y", "y", "Yes", "yes", "YES")) {
        return("No /update performed.")
    }
    
    .update(
        projectName = projectName,
        revisions = revisions,
        token = token,
        ...)
}

.summarize_model_values <- function(model_revs, modelName, model_ids) {
    
    ### Report how many records will be updated
    num_recs <- num_current_recs <- length(model_revs)
    rec_names <- rec_names_current <- names(model_revs)
    
    # Check if any would be be new
    if (! all(rec_names %in% model_ids)) {
        
        rec_names_current <- rec_names[rec_names %in% model_ids]
        
        rec_names_new <- rec_names[!rec_names %in% model_ids]
        num_new <- length(rec_names_new)
        
        ### Summarize for NEW records
        cat("For model \"", modelName, "\", this update() will create ", num_new, " NEW records:\n    ",
            paste0(rec_names_new, collapse = "\n    "),
            "\nWARNING: Check the above carefully. Once created, there is currently no way to remove records from magma.\n",
            sep="")
        
        num_current_recs <- num_recs - num_new
    }
    
    ### Summarize for current records.
    cat("For model \"", modelName, "\", this update() will update ", num_current_recs, " records",
        if (num_current_recs==0) {
                ".\n"
            } else {
                paste0(":\n    ", paste0(rec_names_current, collapse = "\n    "))
            },
        "\n",
        sep="")
    
    model_revs
}


#' Analogous to the '/update' function of magma
#' @inheritParams retrieve
#' @param projectName Single string. The name of the project to upload data to.
#' @param modelName Single string. The name of the model to upload data to.
#' @param attributeName String naming the matrix attribute for which to upload data.
#' @param matrix A matrix or dataframe containing the data to upload to magma.
#' colnames must be record identifiers, and rownames should match the values of 'validation' associated with the target 'attribute'.
#' 
#' Alternatively, the location of a file containing such a data.
#' @param separator String indicating the field separator to use if providing \code{matrix} as a file location.
#' Default = ","
#' @param auto.proceed Logical. When set to TRUE, the function does not ask beforee proceeding forwards to run '/update'.
#' @param ... Additional parameters passed along to the internal `.update()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \itemize{
#' \item \code{request.only} (Logical) & \code{json.params.only} (Logical) which stop the function short and return the values that would have been sent to magma
#' \item \code{verbose} (Logical) sets whether to report the status of the '/update' request to magma.
#' \item \code{url.base} (String) used to direct towards production versus staging versus development versions of magma.
#' }
#' 
#' @return None directly. The function sends data to magma and the only outputs are diagnostic messages.
#' 
#' @details This function utilizes the magma/query function, documented here \url{https://mountetna.github.io/magma.html#update},
#' to upload data to a matrix attribute (named \code{attributeName}) of the \code{modelName} model of \code{projectName} project.
#' 
#' \code{Matrix} data are provided either as a matrix, dataframe, or file path which points toward such data.
#' If given as a file path, the \code{separator} can be used to adjust for whether the file is a csv, \code{separator = ","}, or tsv, \code{separator = "\t"}.
#' 
#' Data is then validated by ensuring that all row names are among the options in the attribute's 'validation', and rows are reordered to be in the same order as these options.
#' 
#' For any missing 'validation' options, NAs are added.
#' 
#' Column names are then checked against record identifiers of the target model.
#' The numbers of new and old record names which will be targeted with the requested update are reported in ordeer to give the user a chance to double-check that the update should proceed.
#' 
#' NOTE: Always check carefully before proceeding. Data can be overwritten with NAs or zeros or the like, but improperly named records cannot be easily removed.
#' 
#' Finally, the data is transformed into the structure required for /update to be called via a curl request,
#' and that request is performed via a call to \code{\link{updateValues}}
#' 
#' @seealso
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying magma/update function.
#' 
#' \code{\link{updateValues}} for a the more general version of this function which can handle non-matrix data types.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     
#'     ### Note that you likely do not have write-permissions for the 'example'
#'     #   project, so this code can be expected to give an authorization error.
#'     
#'     ### Retrieve some data from magma, then update that same data.
#'     mat <- retrieveMatrix("example", "rna_seq", "all", "gene_tpm")
#' 
#'     updateMatrix(
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         attributeName = "gene_tpm",
#'         matrix = mat)
#' }
#'
#' @importFrom utils read.csv
updateMatrix <- function(
    projectName,
    modelName,
    attributeName,
    matrix,
    separator = ",",
    auto.proceed = FALSE,
    token = .get_TOKEN(),
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
    temp <- retrieveTemplate(projectName, token = token, ...)
    row_options <- 
        temp$models[[modelName]]$template$attributes[[attributeName]]$validation$value
    
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
    
    ### Convert & Pass to updateValues() to:
    # Summarize
    # Check
    # Upload
    
    # Transform data into a nested list
    # Note: Do not supply recordNames directly to vapply as any "-" will be
    #   converted to "."
    rec_att_vals <- vapply(seq_len(ncol(matrix)), function(x) {
        x <- list(matrix[,x])
        names(x) <- attributeName
        list(x)
        }, FUN.VALUE = list(1))
    names(rec_att_vals) <- colnames(matrix)
    
    revs <- list(rec_att_vals)
    names(revs) <- modelName
    
    ### Upload
    updateValues(
        projectName = projectName,
        revisions = revs,
        token = token,
        auto.proceed = auto.proceed,
        ...)
}

.ask_before_proceeding <- function(auto.proceed = FALSE) {
    
    if (!auto.proceed) {
        
        if (!interactive()) {
            warning("To run in non-interactive mode, set 'auto.proceed = TRUE'.")
        }
        
        return(readline("Proceed, Y/n? "))
    } else {
        return("Y")
    }
}

.update <- function(
    projectName,
    revisions,
    request.only = FALSE,
    json.params.only = FALSE,
    url.base = .get_URL(),
    token = .get_TOKEN(),
    verbose = TRUE
) {
    
    ### Put together the request 
    jsonParams <- list(
        project_name = projectName,
        revisions = revisions)
    
    requestBody <- jsonlite::toJSON(jsonParams, auto_unbox = TRUE)
    
    ### Output here if requested.
    if (request.only) {
        return(requestBody)
    }
    if (json.params.only) {
        return(jsonParams)
    }
    
    ### Perform '\update'
    curl_out <- .perform_curl_get(
        fxn = "/update",
        requestBody, token, url.base,
        parse = FALSE, verbose = verbose)
}
