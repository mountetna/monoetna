#' Analagous to the '/query' fucntion of magma
#' @inheritParams retrieve
#' @param queryTerms An array of query predicates, see \url{https://mountetna.github.io/magma.html#query}
#' @param format Either "list" or "df" (=dataframe), this sets the desired output format
#' @param ... Additional parameters passed along to the internal `.query()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{request.only} (Logical), \code{json.params.only} (Logical), \code{verbose} (Logical), or \code{url.base} (String which can be user to direct toward production versus staging versus development versions of maagma).
#' @return Either list conversion of the raw '/query' JSON output if \code{format == "list"}, default,
#' 
#' OR a dataframe conversion if \code{format = "df"}
#' @export
#' @examples
#' 
#' ##### FILLER
#' 
#' # Unless a working TOKEN is hard-coded, or it is in an interactive session,
#' #   this code will not work.
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token.
#'     retrieveJSON(
#'         projectName = "ipi",
#'         modelName = "patient",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
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
    url.base = "https://magma.ucsf.edu",
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
    
    h <- RCurl::basicTextGatherer()
    
    RCurl::curlPerform(
        url = paste0(url.base,"/query"),
        httpheader = c('Content-Type' = "application/json", 'Authorization' = paste0('Etna ', token)),
        postfields = requestBody,
        writefunction = h$update,
        verbose = verbose
        )
    
    if (h$value() == "You are unauthorized") {
        stop("Unauthorized. Are you signed into vpn? If yes, run `rm(.TOKEN)` then retry.")
    }
    
    ### Output
    jsonlite::fromJSON(h$value())
}
