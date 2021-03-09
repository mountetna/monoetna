#' Download data from magma that is stored as a matrix
#' @inheritParams retrieve
#' @return a matrix
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     retrieveMatrix(
#'         target = magma,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         recordNames = "all",
#'         attributeNames = "gene_counts")
#' }
#' 
retrieveMatrix <- function(
    target,
    projectName,
    modelName,
    recordNames = "all",
    attributeNames,
    filter = "",
    page = NULL,
    pageSize = 10,
    ...
) {
    
    # Break recordNames into chunks of 10 or fewer
    if (identical(recordNames, "all")) {
        recordNames <- retrieveIds(
            target, projectName, modelName)
    }
    
    sets <- split(recordNames, ceiling(seq_along(recordNames)/10))
    
    # Pull data for each chunk individually, and collate
    chunk_data <- lapply(sets,
            function(x) {
                new <- .matrix_retrieval_chunk(
                    target,
                    projectName = projectName,
                    modelName = modelName,
                    recordNames = x,
                    attributeNames = attributeNames,
                    filter = filter,
                    page = page,
                    pageSize = pageSize,
                    ...)
           })
    data <- do.call(cbind, chunk_data)
    
    # Add row names
    template <- retrieveTemplate(target, projectName = projectName)
    rownames <- template$models[[modelName]]$template$attributes[[attributeNames]]$options
    
    rownames(data) <- rownames
    
    data
}

.matrix_retrieval_chunk <- function(
    target,
    projectName,
    modelName,
    recordNames,
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    request.only = FALSE,
    ...
) {
    
    if (!length(attributeNames) == 1 || attributeNames == "all") {
        stop("This function only works for one attribute at a time.")
    }
    
    # Retrieve as json to get the matrix data
    json <- retrieveJSON(
        target = target,
        projectName = projectName,
        modelName = modelName,
        recordNames = recordNames,
        attributeNames = attributeNames,
        filter = filter,
        page = page,
        pageSize = pageSize,
        request.only = request.only,
        hideTemplate = TRUE,
        ...)
    
    # Extract matrix data as a list of columns
    data_cols <- lapply(
        recordNames,
        function(x) {
            (json$models[[modelName]]$documents[[x]])[[attributeNames]]
        })
    
    # Identify any empty records
    empty <- vapply(data_cols, FUN = length, FUN.VALUE = integer(1)) == 0
    if (any(empty)) {
        for (record in recordNames[empty]){
            warning("Empty record, ", record,", was ignored.")
        }
    }
    
    # Convert from list to matrix
    data <- as.matrix(do.call(data.frame, data_cols[!empty]))
    
    # Add column names
    colnames(data) <- recordNames[!empty]
    
    data
}
