#' Helper functions & special cases of magma /retrieve
#' @name retrieve_SpecialCases
#' @inheritParams retrieve
#' @details
#' These functions aim to help users determine acceptable inputs to other magmaR functions without needing to leave R.
#' 
#' They make properly crafted calls to \code{\link{retrieve}} which target either the "template" or "identifier" special cases outlined in \url{https://mountetna.github.io/magma.html#retrieve},
#' followed by directly returning the output (\code{retrieveTemplate} and \code{retrieveIds}),
#' by returning just a targeted portion of that output (\code{retrieveModels}),
#' or by returning a targeted portion of a subsequent single-record call to \code{\link{retrieve}} (\code{retrieveAttributes}).
#' @return
#' retrieveTemplate = a list conversion of the project's template json.
#' 
#' retrieveModels = a string vector of model names
#' 
#' retrieveIds = a string vector of record names/identifiers.
#' 
#' retrieveAttributes = a string vector of attribute names.
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     
#'     retrieveTemplate(
#'         projectName = "example")
#'         
#'     retrieveModels(
#'         projectName = "example")
#'         
#'     retrieveIds(
#'         projectName = "example",
#'         modelName = "rna_seq")
#'         
#'     retrieveAttributes(
#'         projectName = "example",
#'         modelName = "subject")
#' }
NULL

#' @describeIn retrieve_SpecialCases Retrieve the template for a given project
#' @export
retrieveTemplate <- function(
    target,
    projectName,
    ...
) {
    
    .retrieve(
        target = target,
        projectName = projectName,
        modelName = "all",
        recordNames = "[]",
        attributeNames = "all",
        format = "json",
        ...)
}

#' @describeIn retrieve_SpecialCases Retrieve the modelNames for a given project
#' @export
retrieveModels <- function(
    tagret,
    projectName,
    ...
) {
    names(retrieveTemplate(target, projectName, ...)$models)
}

#' @describeIn retrieve_SpecialCases Retrieve all the identifiers/recordNames for a given project-model pair.
#' @export
retrieveIds <- function(
    target,
    projectName,
    modelName,
    ...
) {
    
    ### This is very inefficient and should be fixed.
    # This chunk is all that should be required
    .retrieve(
        target = target,
        projectName = projectName,
        modelName = modelName,
        recordNames = "all",
        attributeNames = "identifier",
        format = "tsv",
        ...)[,,drop = TRUE]
}

#' @describeIn retrieve_SpecialCases Retrieve all the attribute options for a given project-model pair.
#' @export
retrieveAttributes <- function(
    target,
    projectName,
    modelName,
    ...
) {
    
    rec <- retrieveIds(target, projectName, modelName, ...)[1]

    .retrieve(
        target, projectName, modelName, recordNames = rec,
        attributeNames = "all", format = "tsv",
        names.only = TRUE, ...)
}
