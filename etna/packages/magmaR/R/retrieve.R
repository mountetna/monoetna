#' Download data from magma as a tsv, and convert to a data.frame
#' @description Analogous to the '/retrieve' function of magma, with \code{format = "tsv"}
#' @param target A list, which can be created using \code{\link{magmaRset}}, containing your authorization 'token' (a string), a 'url' of magma to target (a string), and optional 'opts' for specifying additions parameters for curl requests (a named list).
#' @param projectName Single string. The name of the project you would like to interact with. For options, see \code{\link{retrieveProjects}}.
#' @param modelName Single string. The name of the subset data structure within the project, which are referred to as 'model's in magma, to interact with.
#' For options, see \code{\link{retrieveModels}} or https://timur.ucsf.edu/<projectName>/map.
#' @param recordNames Single string or string vector indicating which particular sample/tube/etc. records to target.
#' Options are "all" or any combination of individual record names. To retrieve individual options, see \code{\link{retrieveIds}}.
#' @param attributeNames Single string or string vector indicating which features of the data to target.
#' Options are "all" or any combination of individual attribute names. To retrieve individual options, see \code{\link{retrieveAttributes}}.
#' @param filter String. Potential filter(s) of the data.
#' Example: "<targetAttributeName>~GYN" to filter to records where <targetAttributeName> contains "GYN".
#' 
#' Refer to \url{https://mountetna.github.io/magma.html#retrieve} for more details about options and format.
#' @param showDisconnected Boolean. Set to true to access "disconnected" records, which are records that are missing an (upstream) parent linkage, and so do not connect up with the project's top-level project record.
#' Generally, disconnected records are ones that were deemed low quality in some way, thus were purposefully disconnected from the rest of the dataset.
#' But sometimes a record might just be disconnected because an upload went awry.
#' @param pageSize Integer. For retrieving just a portion of the data, sets slice/page size, which is equivalent to the a number of rows.
#' @param page Integer. For retrieving just a portion of the data, sets which slice to get.
#' @param ... Additional parameters passed along to the internal `.retrieve()`, `.query()`, or `.update()` functions,
#' for troubleshooting or advanced-user purposes only: \itemize{
#' \item \code{request.only} (Logical) & \code{json.params.only} (Logical) which 1) stop the function before its main curl request to magma and 2) returns the values that would have been sent to magma in either of two formats.
#' \item \code{verbose} (Logical) sets whether to report the status of the curl request after it is performed.
#' }
#' @return A dataframe
#' @details This function makes a curl get request to magma/retrieve, with properly reformatted versions of user inputs, plus \code{format = "tsv"}.
#' Then, it converts the tsv-string output into a dataframe.
#' 
#' Note: When \code{format = "tsv"}, magma/retrieve returns just an identifier for matrix-type attributes.
#' To retrieve underlying data for such attributes, use the specialized \code{\link{retrieveMatrix}} function.
#' @seealso
#' \code{\link{retrieveMatrix}} for retrieving attributes of type matrix.
#' 
#' \code{\link{retrieveJSON}} for similar functionality to \code{retrieve}, but where the call to magma/retrieve is made with \code{format = "json"} and the output is a list.
#' This output often contains more information, and can retrieve data for attribute types of type matrix, which are not returned by the current function.
#' But in most cases, the data returned by \code{retrieve} and \code{retrieveMatrix} will suffice.
#' 
#' \url{https://mountetna.github.io/magma.html#retrieve} for documentation of the underlying magma/retrieve function.
#' 
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     # Now we can retrieve data with...
#'     retrieve(
#'         target = magma,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
#' }
#' 
retrieve <- function(
    target,
    projectName,
    modelName,
    recordNames = "all",
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    showDisconnected = FALSE,
    ...
) {
    .retrieve(
        target = target,
        projectName = projectName,
        modelName = modelName,
        recordNames = recordNames,
        attributeNames = attributeNames,
        filter = filter,
        format = "tsv",
        page = page,
        pageSize = pageSize,
        showDisconnected = showDisconnected,
        ...)
}

#' Download data from magma as a json, and convert to a list
#' @description Analogous to the '/retrieve' function of magma, with \code{format = "json"}
#' @inheritParams retrieve
#' @param hideTemplate Logical. Allows to leave out the project template from the return. Often this does not matter much, but the template can be bulky.
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
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     # Now we can retrieve data as json (->list) with...
#'     json_out <- retrieveJSON(
#'         target = magma,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
#'     # The return will be a nested list with data in a 'documents' element and
#'     #  some extra information about each attribute in a 'template' element.
#'     str(json_out, max.level = 4)
#'
#'     # Often, the 'template' part is bulky but not needed, so its retrieval may
#'     #  be turned off by giving hideTemplaate = TRUE'
#'     json_out <- retrieveJSON(
#'         target = magma,
#'         projectName = "example",
#'         modelName = "rna_seq",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "",
#'         hideTemplate = TRUE)
#'     str(json_out, max.level = 4)
#' }
#' 
retrieveJSON <- function(
    target,
    projectName,
    modelName,
    recordNames = "all",
    attributeNames = "all",
    filter = "",
    page = NULL,
    pageSize = 10,
    showDisconnected = FALSE,
    hideTemplate = FALSE,
    ...
) {
    
    .retrieve(
        target = target,
        projectName = projectName,
        modelName = modelName,
        recordNames = recordNames,
        attributeNames = attributeNames,
        filter = filter,
        format = "json",
        page = page,
        pageSize = pageSize,
        showDisconnected = showDisconnected,
        hideTemplate = hideTemplate,
        ...)
}

.retrieve <- function(
    target,
    projectName,
    modelName,
    recordNames,
    attributeNames,
    format = c("tsv", "json"),
    filter = "",
    page = NULL,
    pageSize = 10,
    showDisconnected = FALSE,
    hideTemplate = FALSE,
    names.only = FALSE,
    request.only = FALSE,
    json.params.only = FALSE,
    raw.return = FALSE,
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
        show_disconnected = showDisconnected,
        hide_templates = hideTemplate,
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
        target, requestBody,
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
