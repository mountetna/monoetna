#' A matrix-specific wrapper of \code{\link{updateValues}}
#' @description A matrix-specific wrapper of \code{\link{updateValues}} which can take in a matrix, data.frame, or file path, directly.
#' @inheritParams retrieve
#' @param attributeName String naming the matrix attribute for which to upload data.
#' @param matrix A matrix or dataframe containing the data to upload to magma.
#' 
#' Alternatively, a String specifying the file path of a file containing such data.
#' 
#' No matter the provision method, colnames must be record identifiers, and rownames should match the values of 'options' associated with the target 'attribute'.
#' Check the 'See Also' section below for how to determine the needed 'options'.
#' @param separator String indicating the field separator to use if providing \code{matrix} as a file path.
#' Default = \code{","}.
#' @param auto.proceed Logical. When set to TRUE, the function does not ask before proceeding forward with the 'magma/update'.
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
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying magma/update function.
#' 
#' \code{\link{updateValues}} for a the more general version of this function which can handle non-matrix data types.
#' 
#' \code{\link{retrieveTemplate}}, then check the \code{<output>$models$<modelName>$template$attributes$<attributeName>$options} to explore the rownames that your matrix should have.
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
    separator = ",",
    auto.proceed = FALSE,
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
    temp <- retrieveTemplate(target, projectName)
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
        target = target,
        projectName = projectName,
        revisions = revs,
        auto.proceed = auto.proceed,
        ...)
}
