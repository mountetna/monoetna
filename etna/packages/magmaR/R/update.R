#' Analogous to the '/update' function of magma
#' @inheritParams retrieve
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
#' @param ... Additional parameters passed along to the internal `.update()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{request.only} (Logical), \code{json.params.only} (Logical), \code{verbose} (Logical), or \code{url.base} (String which can be used to direct toward production versus staging versus development versions of magma).
#' @return None directly.  The function sends data to magma, and the only output is whether that send was successful.
#' @details This function initially mimics the activity of the magma/query which is documented here \url{https://mountetna.github.io/magma.html#query}.
#' 
#' Afterwards, the json list output of magma/query is converted into an R list, and then the \code{format} input determines whether it should be wrangled further:
#' \itemize{
#' \item \code{format = "list"}, default: R list output directly.
#' \item \code{format = "df"}: R list converted into a dataframe where data comes from the list$answer and column names come from the list$format
#' }
#' @seealso
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying magma/update function.
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

updateMagma <- function(
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
