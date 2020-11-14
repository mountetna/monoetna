#' Download data from magma as a tsv and convert to a data.frame
#' @param token your TOKEN from https://janus.ucsf.edu
#' @param projectName the name of the project you would like to download data from, exactly as it appears at https://timur.ucsf.edu
#' @param modelName the name of the piece of the magma model form which to download data
#' @param recordNames,attributeNames which records and attributes to grab
#' @param filter potential filter of the data
#' @param pageSize,page Integers. For retrieving just a portion of the data, sets slice/page size which is equivalent to a number of rows, and which slice to get.
#' @param connected.only Logical. Whether data without a parent should be retained. 
#' @param ... Additional parameters passed along to the internal `.retrieve()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{request.only} (Logical), \code{json.params.only} (Logical), \code{verbose} (Logical), or \code{url.base} (String which can be user to direct toward production versus staging versus development versions of maagma).
#' @return A dataframe representation of the requested data (or a list if the data was not rectangular).
#' @export
#' @examples
#' 
#' # Unless a working TOKEN is hard-coded, or it is in an interactive session,
#' #   this code will not work.
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token.
#'     retrieve(
#'         projectName = "ipi",
#'         modelName = "patient",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
#' }
#' 
#' @importFrom utils URLencode

retrieve <- function(
    projectName,
    modelName,
    recordNames = "all",
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    connected.only = TRUE,
    token = .get_TOKEN(),
    ...
) {
    .retrieve(
        projectName = projectName,
        modelName = modelName,
        recordNames = recordNames,
        attributeNames = attributeNames,
        filter = filter,
        format = "tsv",
        page = page,
        pageSize = pageSize,
        connected.only = TRUE,
        token = token,
        ...)
}

#' Download data from magma as a json, output as a json string (list)
#' @inheritParams retrieve
#' @param filter potential filter of the data
#' @return a json string or list of such strings
#' @export
#' @examples
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
retrieveJSON <- function(
    projectName,
    modelName,
    recordNames = "all",
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    token = .get_TOKEN(),
    ...
) {
    
    .retrieve(
        projectName = projectName,
        modelName = modelName,
        recordNames = recordNames,
        attributeNames = attributeNames,
        filter = filter,
        format = "json",
        page = page,
        pageSize = pageSize,
        token = token,
        ...)
}

.retrieve <- function(
    projectName,
    modelName,
    recordNames,
    attributeNames,
    format = c("tsv", "json"),
    filter = "",
    page = NULL,
    pageSize = 10,
    connected.only = TRUE,
    names.only = FALSE,
    request.only = FALSE,
    json.params.only = FALSE,
    raw.return = FALSE,
    url.base = "https://magma.ucsf.edu",
    token = .get_TOKEN(),
    verbose = FALSE
) {
    
    format <- match.arg(format)
    
    ### Put together the request 
    
    jsonParams <- list(
        project_name = projectName,
        model_name = modelName,
        record_names = .as_array_unless_all_or_brackets(recordNames),
        attribute_names = .as_array_unless_all_or_identifier(attributeNames),
        filter = filter,
        format = format)

    if (!identical(page, NULL)) {
        jsonParams$page <- page
        jsonParams$page_size <- pageSize
    }
    
    requestBody <- jsonlite::toJSON(jsonParams, auto_unbox = TRUE)
    
    ### Output here if requested.
    if (request.only) {
        return(requestBody)
    }
    if (json.params.only) {
        return(jsonParams)
    }
    
    ### Perform '\retrieve'-al
    
    h <- RCurl::basicTextGatherer()
    
    RCurl::curlPerform(
        url = paste0(url.base,"/retrieve"),
        httpheader = c('Content-Type' = "application/json", 'Authorization' = paste0('Etna ', token)),
        postfields = requestBody,
        writefunction = h$update,
        verbose = verbose
        )
    
    if (h$value() == "You are unauthorized") {
        stop("Unauthorized. Are you signed into vpn? If yes, run `rm(.TOKEN)` then retry.")
    }
    
    ### Output
    if (raw.return) {
        h$value()
    } else {
        if (format=="tsv") {
            .parse_tsv(h$value(), names.only, connected.only)
        } else {
            jsonlite::fromJSON(h$value())
        }
    }
}
