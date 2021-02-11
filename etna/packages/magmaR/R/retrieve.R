#' Download data from magma as a tsv, and convert to a data.frame
#' @param token your personal TOKEN from \url{https://janus.ucsf.edu}. 
#' 
#' When not explicitly given in a function call: you will be prompted to input token, one time.
#' This user provided token will then be stored as a hidden variable, \code{.TOKEN}, in the global R environment,
#' and all future magmaR calls without an explicitly provided token will turn to this \code{.TOKEN}.
#' @param projectName Single string. The name of the project you would like to download data from, exactly as it appears at https://timur.ucsf.edu
#' @param modelName Single string. The name of data structure within your project, which are referred to as 'model's in magma, form which to download data.
#' For options, see \code{\link{retrieveModels}} or https://timur.ucsf.edu/<projectName>/map.
#' @param recordNames Single string or string vector. Which particular sample/tube/etc. records to grab data for.
#' Options are "all" or any combination of individual record names. To retrieve individual options, see \code{\link{retrieveIds}}.
#' @param attributeNames Single string or string vector. Which features of the data to obtain.
#' Options are "all" or any combination of individual attribute names. To retrieve individual options, see \code{\link{retrieveAttributes}}.
#' @param filter String. Potential filter of the data.
#' Example: "<targetAttributeName>~GYN" to filter to records where <targetAttributeName> contains "GYN".
#' 
#' Refer to \url{https://mountetna.github.io/magma.html#retrieve} for more details.
#' @param pageSize,page Integers. For retrieving just a portion of the data, sets slice/page size, which is equivalent to the a number of rows, and which slice to get.
#' @param connected.only \emph{Not implemented yet.} Logical. Whether data without a parent should be retained. 
#' @param ... Additional parameters passed along to the internal `.retrieve()` function.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{request.only} (Logical), \code{json.params.only} (Logical), \code{verbose} (Logical), or \code{url.base} (String which can be used to direct toward production versus staging versus development versions of magma).
#' @return A dataframe
#' @details This function makes a call to magma/retrieve with \code{format = "tsv"}.
#' Then, it converts the tsv-string output into a dataframe.
#' @seealso
#' \url{https://mountetna.github.io/magma.html#retrieve} for documentation of the underlying magma/retrieve function.
#' 
#' \code{\link{retrieveModels}}, \code{\link{retrieveIds}}, and \code{\link{retrieveAttributes}} for exploring options for the \code{modelName}, \code{recordNames}, and \code{attributeNames} inputs, respectively.
#' 
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     retrieve(
#'         projectName = "ipi",
#'         modelName = "patient",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
#' }
#' 
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

#' Download data from magma as a json, and view it as a list
#' @inheritParams retrieve
#' @param filter potential filter of the data
#' @return A list
#' @details This function makes a call to magma/retrieve with \code{format = "json"}.
#' Then, it converts the json output into a list which is more compatible with R.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
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
    url.base = .get_URL(),
    token = .get_TOKEN(),
    verbose = FALSE
) {
    
    format <- match.arg(format)
    
    ### Put together the request 
    jsonParams <- list(
        project_name = projectName,
        model_name = modelName,
        record_names = .I_list_unless_all(recordNames),
        attribute_names = .I_list_unless_all_or_identifier(attributeNames),
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
    curl_out <- .perform_curl_get(
        fxn = "/retrieve",
        requestBody, token, url.base,
        parse = TRUE, verbose = verbose)
    
    ### Output
    if (raw.return) {
        curl_out
    } else {
        if (format=="tsv") {
            .parse_tsv(curl_out, names.only, connected.only)
        } else {
            jsonlite::fromJSON(curl_out)
        }
    }
}

.I_list_unless_all <- function(values) {
    
    if (identical(values,"[]")) {
        return(I(list()))
    }
    
    ifelse(test = identical(values, "all"),
           yes = values,
           no = I(as.list(values)))
}

.I_list_unless_all_or_identifier <- function(values) {
    
    ifelse(test = (identical(values, "all") || identical(values, "identifier")),
           yes = values,
           no = I(as.list(values)))
}
