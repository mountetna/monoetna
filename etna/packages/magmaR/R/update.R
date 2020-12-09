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
#' @details This function mimics the activity of the magma/query function, documented here \url{https://mountetna.github.io/magma.html#update},
#' with the main difference being that the \code{revisions} input should be in nested-list format.
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
#'         projectName = projectName,
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
    ...
) {
    
    .update(
        projectName = projectName,
        revisions = revisions,
        token = token,
        ...)
}

#' Analogous to the '/update' function of magma
#' @inheritParams retrieve
#' @param projectName Single string. The name of the project to upload data to.
#' @param modelName Single string. The name of the model to upload data to.
#' @param attributeName String naming the matrix attribute for which to upload data.
#' @param matrix A matrix or dataframe containing the data to upload to magma.
#' colnames must be record identifiers, and rownames should match the 'options' associated with the target 'attribute'.
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
#' Data is then validated by ensuring that all row names are among the attribute's 'options', and rows are reordered to be in thee same order as these 'options'.
#' Column names are then checked against record identifiers of the target model.
#' The numbers of new and old record names which will be targeted with the requested update are reported in ordeer to give the user a chance to double-check that the update should proceed.
#' 
#' NOTE: Always check carefully before proceeding. Data can be overwritten with NAs or zeros, but improperly named records cannot be easily removed.
#' 
#' Finally, the data is transformed into the structure required for /update to be called via a curl request.
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
        matrix <- read.csv(matrix, sep = separator)
    }
    if (is.null(colnames(matrix))) {
        stop("Colnames of matrix should be record names.")
    }
    
    ### Validate rownames (genes / value names)
    # Obtain rownames
    temp <- retrieveTemplate(projectName, token = token, ...)
    row_options <- 
        temp$models[[modelName]]$template$attributes[[attributeName]]$options
    
    if (length(row_options)<1) {
        warning("Attribute does not have 'options' info: no validation can be performed.")
    } else {
        # Check rownames against options
        if (! all(rownames(matrix) %in% row_options)) {
            stop("Validation error: rownames of 'matrix' are not valid 'options' for ", attributeName,".")
        } else {
            # Add NAs for any 'options' not in the current matrix
            for (not_there in setdiff(row_options, rownames(matrix))) {
                matrix <- rbind(matrix, NA)
                rownames(matrix)[nrow(matrix)] <- not_there
            }
            # Reorder as in options
            matrix <- matrix[row_options, ,drop = FALSE]
        }
    }
    
    ### Report what will happen
    num_updates <- ncol(matrix)
    
    # Check if colnames (record identifiers) would be new
    ids <- retrieveIds(projectName, modelName, token, ...)
    if (! all(colnames(matrix) %in% ids)) {
        new <- colnames(matrix)[!colnames(matrix) %in% ids]
        num_new <- length(new)
        
        cat("This update() will create", num_new, "new records:",
            paste0(head(new,3), collapse = ", "),"...\n")
        
        num_updates <- num_updates - num_new
    }
    
    cat("This update() will update data for", num_updates, "records.\n")
    
    ### Check if should move forward
    .ask_before_proceeding(auto.proceed)
    
    ### Transform data into a nested list
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
    .update(
        projectName = projectName,
        revisions = revs,
        token = token,
        ...)
}

#' @importFrom utils head
.ask_before_proceeding <- function(auto.proceed = FALSE) {
    
    if (!auto.proceed) {
        
        if (!interactive()) {
            stop("To run in non-interactive mode, set 'auto.proceed = TRUE'.")
        }
        
        go <- readline("Proceed, Y/n? ")
        if (! go %in% c("Y", "y", "Yes", "yes", "YES")) {
            stop("User-requested exit.", call. = FALSE)
        }
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
