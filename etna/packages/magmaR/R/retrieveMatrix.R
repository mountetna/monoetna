#' Download data stored as a matrix
#' @inheritParams retrieve
#' @return a matrix
#' @export
#' @examples
#' 
#' # Unless a working TOKEN is hard-coded, or it is in an interactive session,
#' #   this code will not work.
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token on first.
#'     
#'     ids <- retrieve("ipi", "rna_seq", attributeNames = "identifier")
#'     
#'     retrieveMatrix(
#'         projectName = "ipi",
#'         modelName = "rna_seq",
#'         recordNames = ids[1:4],
#'         attributeNames = "gene_counts")
#' }
#' 
retrieveMatrix <- function(
    projectName,
    modelName,
    recordNames,
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    token = .get_TOKEN(),
    ...
) {
    
    # Break recordNames into chunks of 10 or fewer
    if (identical(recordNames, "all")) {
        recordNames <- .retrieve(
            projectName, modelName, attributeNames = "identifier", ...)
    }
    
    sets <- split(recordNames, ceiling(seq_along(recordNames)/10))
    
    # Pull data for each chunk individually, and collate
    chunk_data <- lapply(sets,
           function(x) {
               new <- .matrix_retrieval_chunk(
                    projectName = projectName,
                    modelName = modelName,
                    recordNames = x,
                    attributeNames = attributeNames,
                    filter = filter,
                    page = page,
                    pageSize = pageSize,
                    token = token,
                    ...)
           })
    data <- do.call(cbind, chunk_data)
    
    # Add row names
    template <- retrieveTemplate(projectName = projectName, token = token)
    rownames <- template$models[[modelName]]$template$attributes[[attributeNames]]$options
    
    rownames(data) <- rownames
    
    data
}

.matrix_retrieval_chunk <- function(
    projectName,
    modelName,
    recordNames,
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    request.only = FALSE,
    token = .get_TOKEN(),
    ...
) {
    
    if (!length(attributeNames) == 1 || attributeNames %in% c("all", "identifier")) {
        stop("This function only works for one attribute at a time.")
    }
    
    json <- retrieveJSON(
        projectName = projectName,
        modelName = modelName,
        recordNames = recordNames,
        attributeNames = attributeNames,
        filter = filter,
        page = page,
        pageSize = pageSize,
        request.only = request.only,
        token = token,
        ...)
    
    data <- sapply(
        recordNames,
        function(x) {
            (json$models[[modelName]]$documents[[x]])[[attributeNames]]
        })
    colnames(data) <- recordNames
    
    data
}
