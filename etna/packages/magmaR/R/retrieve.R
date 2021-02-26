#' Download data from magma as a tsv, and convert to a data.frame
#' @param token Single string. Your personal token from \url{https://janus.ucsf.edu}. 
#' 
#' When not explicitly given in a function call: you will be prompted to input your token, one time.
#' This user provided token will then be stored as a hidden variable, \code{.MAGMAR_TOKEN}, in the global R environment,
#' and all future magmaR calls without a \code{token} explicitly provided will turn to this \code{.MAGMAR_TOKEN}.
#' @param projectName Single string. The name of the project you would like to interact with. For options, see \code{\link{retrieveProjects}}.
#' @param modelName Single string. The name of the subset data structure within the project, which are referred to as 'model's in magma, to interact with.
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
#' @param ... Additional parameters passed along to internal `.retrieve()`, `.query()`, or `.update()` functions,
#' for troubleshooting or advanced-user purposes only: \itemize{
#' \item \code{request.only} (Logical) & \code{json.params.only} (Logical) which stop the function before it performs any curl requests and instead outputs the values that would have been sent to magma in, either of two formats.
#' \item \code{verbose} (Logical) sets whether to report the status of the curl request after it is performed.
#' \item \code{url.base} (String) used to direct towards production versus staging versus development versions of magma. See \code{\link{magma-environments}}
#' }
#' @return A dataframe
#' @details This function makes a curl get request to magma/retrieve, with properly reformatted versions of user inputs, plus \code{format = "tsv"}.
#' Then, it converts the tsv-string output into a dataframe.
#' 
#' Note: When \code{format = "tsv"}, magma/retrieve returns just an identifier for matrix-type attributes.
#' To retrieve underlying data for such attributes, use \code{\link{retrieveMatrix}} which is a specialized wrapper around \code{\link{retrieveJSON}}.
#' @seealso
#' \code{\link{retrieveMatrix}} for retrieving matrix data types.
#' 
#' \code{\link{retrieveJSON}} for similar functionality, but where the call to magma/retrieve is made with \code{format = "json"} and the output is a list.
#' 
#' \url{https://mountetna.github.io/magma.html#retrieve} for documentation of the underlying magma/retrieve function.
#' 
#' \code{\link{retrieveProjects}} for exploring options for the \code{projectName} input.
#' 
#' \code{\link{retrieveModels}}, \code{\link{retrieveIds}}, and \code{\link{retrieveAttributes}} for exploring options for the \code{modelName}, \code{recordNames}, and \code{attributeNames} inputs, respectively.
#' 
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     retrieve(
#'         projectName = "example",
#'         modelName = "rna_seq",
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
        token = token,
        ...)
}

#' Download data from magma as a json, and convert to a list
#' @inheritParams retrieve
#' @return A list
#' @details This function makes a call to magma/retrieve with \code{format = "json"}.
#' Then, it converts the json output into a list which is more compatible with R.
#' @seealso
#' \code{\link{retrieve}} for similar functionality, but where the call to magma/retrieve will be made with \code{format = "tsv"} and the output is a dataframe.
#' 
#' \code{\link{retrieveMatrix}} for matrix data-targeted utilization of this current \code{retreiveJSON} function, followed by automated restructuring of the return into a matrix format.
#' 
#' \url{https://mountetna.github.io/magma.html#retrieve} for documentation of the underlying magma/retrieve function.
#' 
#' \code{\link{retrieveProjects}} for exploring options for the \code{projectName} input.
#' 
#' \code{\link{retrieveModels}}, \code{\link{retrieveIds}}, and \code{\link{retrieveAttributes}} for exploring options for the \code{modelName}, \code{recordNames}, and \code{attributeNames} inputs, respectively.
#' 
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     retrieveJSON(
#'         projectName = "example",
#'         modelName = "rna_seq",
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
        record_names = .match_expected_recName_structure(recordNames),
        attribute_names = .match_expected_attName_structure(attributeNames),
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
            .parse_tsv(curl_out, names.only)
        } else {
            jsonlite::fromJSON(curl_out)
        }
    }
}

.match_expected_recName_structure <- function(values) {
    
    if (identical(values,"[]")) {
        return(I(list()))
    }
    
    .match_expected_common(values)
}

.match_expected_attName_structure <- function(values) {
    
    if (identical(values, "identifier")) {
        return(values)
    }
    
    .match_expected_common(values)
}

.match_expected_common <- function(values) {
    # all -> all
    # vector -> I(vector)
    # single -> I(list(single))
    
    if (identical(values, "all")) {
        return(values)
    }
    
    if (length(values) > 1) {
        I(values)
    } else {
        I(list(values))
    }
}
