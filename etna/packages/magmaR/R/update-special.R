#' A matrix-specific wrapper of updateValues()
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
#' @param auto.proceed Logical. When set to TRUE, the function does not ask before proceeding forward to run '/update'.
#' @param ... Additional parameters passed along to the internal `.update()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \itemize{
#' \item \code{request.only} (Logical) & \code{json.params.only} (Logical) which stop the function short and return the values that would have been sent to magma
#' \item \code{verbose} (Logical) sets whether to report the status of the '/update' request to magma.
#' \item \code{url.base} (String) used to direct towards production versus staging versus development versions of magma.
#' }
#' 
#' @return None directly.
#' The function sends data to magma, and the only output is whether that send was successful, when \code{verbose = TRUE},
#' or the string "No /update performed." if the user chooses not to proceed with performing the update.
#' @details This function utilizes the magma/query function, documented here \url{https://mountetna.github.io/magma.html#update},
#' to upload data to a matrix attribute (named \code{attributeName}) of the \code{modelName} model of \code{projectName} project.
#' 
#' \code{matrix} data are provided either as a matrix, dataframe, or file path which points toward such data.
#' If given as a file path, the \code{separator} can be used to adjust for whether the file is a csv, \code{default, separator = ","}, or tsv, \code{separator = "\t"}.
#' 
#' Data is then validated by ensuring that all row names are among the options in the attribute's 'validation', and rows are reordered to be in the same order as these options.
#' 
#' For any missing 'validation' options, NAs are added.
#' 
#' The data is then transformed and passed along to \code{\link{updateValues}}.
#' 
#' The \code{updateValues()} function will summarized records to be updated and allow the user to double-check this information before proceeding.
#' 
#' This user-prompt step can be bypassed (useful when running in a non-interactive way) by setting \code{auto.proceed = TRUE}, but NOTE:
#' It is a good idea to always check carefully before proceeding, if possible.
#' Data can be overwritten with NAs or zeros or the like, but improperly named records cannot be easily removed.
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
#'     # project, so this code can be expected to give an authorization error.
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
    names(revs) <- modelName
    
    ### Pass to updateValues() to:
    # Summarize to user
    # Check with the user before proceeding
    # Perform upload
    updateValues(
        projectName = projectName,
        revisions = revs,
        token = token,
        auto.proceed = auto.proceed,
        ...)
}
