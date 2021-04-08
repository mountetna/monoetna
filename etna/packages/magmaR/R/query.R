#' Search-like function that can obtain linked data from distinct models.
#' @description Analogous to the '/query' function of magma.
#' @inheritParams retrieve
#' @param queryTerms A list of strings where list elements are query predicates and verbs. See \url{https://mountetna.github.io/magma.html#query} for details.
#' @param format Either "list" or "df" (=dataframe). This sets the desired output format. The list option is the more raw form.
#' @return A list, default, if \code{format == "list"},
#' 
#' OR A dataframe conversion if \code{format = "df"}
#' @details This function initially mimics the activity of the magma's /query functionality,
#' which is documented here \url{https://mountetna.github.io/magma.html#query}.
#' 
#' Afterwards, the json list output of magma/query is converted into an R list, and then the \code{format} input determines whether it should be wrangled further:
#' \itemize{
#' \item \code{format = "list"}, default: R list output directly.
#' \item \code{format = "df"}: R list converted into a dataframe where data comes from the list$answer and column names come from the list$format
#' }
#' @seealso
#' \url{https://mountetna.github.io/magma.html#query} for documentation of the underlying magma/query function.
#' 
#' \code{\link{retrieveProjects}} for exploring options for the \code{projectName} input.
#' 
#' \code{\link{retrieveModels}}, \code{\link{retrieveIds}}, and \code{\link{retrieveAttributes}} and \code{\link{retrieveTemplate}} for exploring the project structure and determining \code{queryTerm} options.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     ### To obtain the 'group' attribute, from the subject-model, that are
#'     #   associated with records of the rna_seq-model:
#' 
#'     # "Raw" output of query:
#'     query_list <- query(
#'         target = magma,
#'         projectName = "example",
#'         queryTerms = 
#'             list('rna_seq',
#'                  '::all',
#'                  'biospecimen',
#'                  'subject',
#'                  'group'))
#'     print(query_list)
#'                  
#'     # Or instead re-formatted to a dataframe, which may be easier for
#'     #   downstream applications in R:
#'     query_df <- query(
#'         target = magma,
#'         projectName = "example",
#'         queryTerms = 
#'             list('rna_seq',
#'                  '::all',
#'                  'biospecimen',
#'                  'subject',
#'                  'group'),
#'         format = 'df')
#'     print(query_df)
#' }
#'

query <- function(
    target,
    projectName,
    queryTerms = list(),
    format = c("list", "df"),
    ...
) {
    
    format <- match.arg(format)
    
    raw <- .query(
        target = target,
        projectName = projectName,
        queryTerms = queryTerms,
        ...)
    
    if (format == "list") {
        raw
    } else {
        out <- data.frame(raw$answer)
        colnames(out) <- raw$format
        
        out
    }
}

.query <- function(
    target,
    projectName,
    queryTerms = list(),
    request.only = FALSE,
    json.params.only = FALSE,
    verbose = FALSE
) {
    
    ### Put together the request 
    
    jsonParams <- list(
        project_name = projectName,
        query = queryTerms)

    requestBody <- jsonlite::toJSON(jsonParams, auto_unbox = TRUE)
    
    ### Output here if requested.
    if (request.only) {
        return(requestBody)
    }
    if (json.params.only) {
        return(jsonParams)
    }
    
    ### Perform '\query'
    curl_out <- .perform_curl_get(
        fxn = "/query",
        target, requestBody,
        parse = TRUE, verbose = verbose)
    
    ### Output
    jsonlite::fromJSON(curl_out)
}
