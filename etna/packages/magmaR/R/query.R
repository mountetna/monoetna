#' Search-like function that can obtain linked data from distinct models.
#' @description Analogous to the '/query' function of magma.
#' @inheritParams retrieve
#' @param queryTerms A list of strings where list elements are query predicates and verbs. See \url{https://mountetna.github.io/magma.html#query} for details.
#' @param format Either "list", "tsv" or "df" (=dataframe).
#' This sets the desired output format.
#' The list option is the most raw form, whereas setting to "df" or "tsv" will utilize magma/query's 'format="tsv"' option to have data converted into a tsv string format, before it is output.
#' The tsv string will then read into a data.frame.
#' @return A list, default, if \code{format == "list"},
#' 
#' OR A data.frame if \code{format = "df" or "tsv"}
#' @details This function mimics the activity of the magma's /query functionality,
#' which is documented here \url{https://mountetna.github.io/magma.html#query}.
#' 
#' Afterwards, output formats are adjusted slightly for R compatibility or use-ability, depending on the 'format' requested:
#' \itemize{
#' \item \code{format = "list"}, default: json format returned is output as it's nearest cognate, an R list.
#' \item \code{format = "df" or "tsv"}: tsv string format returned is converted into a data.frame.
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
    format = c("list", "df", "tsv"),
    ...
) {
    
    format <- match.arg(format)
    
    .query(
        target = target,
        projectName = projectName,
        queryTerms = queryTerms,
        format = format,
        ...)
}

.query <- function(
    target,
    projectName,
    queryTerms = list(),
    format = "list",
    request.only = FALSE,
    json.params.only = FALSE,
    verbose = FALSE
) {
    
    ### Put together the request 
    
    jsonParams <- list(
        project_name = projectName,
        query = queryTerms)
    if (format %in% c("df", "tsv")) {
        jsonParams$format <- "tsv"
    }

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
    if (format %in% c("df", "tsv")) {
        .parse_tsv(curl_out)
    } else {
        jsonlite::fromJSON(curl_out)
    }
}
