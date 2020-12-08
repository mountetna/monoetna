#' Analagous to the '/query' fucntion of magma
#' @inheritParams retrieve
#' @param queryTerms A string vector where elements are query predicates. See \url{https://mountetna.github.io/magma.html#query} for details.
#' @param format Either "list" or "df" (=dataframe). This sets the desired output format.
#' @param ... Additional parameters passed along to the internal `.query()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{request.only} (Logical), \code{json.params.only} (Logical), \code{verbose} (Logical), or \code{url.base} (String which can be used to direct toward production versus staging versus development versions of magma).
#' @return A list, default, if \code{format == "list"},
#' 
#' OR A dataframe conversion if \code{format = "df"}
#' @details This function initially mimics the activity of the magma/query which is documented here \url{https://mountetna.github.io/magma.html#query}.
#' 
#' Afterwards, the json list output of magma/query is converted into an R list, and then the \code{format} input determines whether it should be wrangled further:
#' \itemize{
#' \item \code{format = "list"}, default: R list output directly.
#' \item \code{format = "df"}: R list converted into a dataframe where data comes from the list$answer and column names come from the list$format
#' }
#' @seealso
#' \url{https://mountetna.github.io/magma.html#query} for documentation of the underlying magma/query function.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     
#'     ### To obtain the sample-model record-identifiers that are linked to
#'     #   records of the rna_seq-model:
#' 
#'     # "Raw" output of query:
#'     query(
#'         projectName = projectName,
#'         queryTerms = 
#'             list('rna_seq',
#'                  '::all',
#'                  'sample',
#'                  '::identifier'))
#'                  
#'     # Or instead re-formatted to a dataframe, which may be easier for
#'     #   downstream applications in R:
#'     query(
#'         projectName = projectName,
#'         queryTerms = 
#'             list('rna_seq',
#'                  '::all',
#'                  'sample',
#'                  '::identifier'),
#'         format = 'df')
#' }
#'

query <- function(
    projectName,
    queryTerms = list(),
    format = c("list", "df"),
    token = .get_TOKEN(),
    ...
) {
    
    format <- match.arg(format)
    
    raw <- .query(
        projectName = projectName,
        queryTerms = queryTerms,
        token = token,
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
    projectName,
    queryTerms = list(),
    request.only = FALSE,
    json.params.only = FALSE,
    url.base = .get_URL(),
    token = .get_TOKEN(),
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
    curl_out <- .perform_curl(
        fxn = "/query", requestBody, token, url.base, verbose)
    
    ### Output
    jsonlite::fromJSON(curl_out)
}
