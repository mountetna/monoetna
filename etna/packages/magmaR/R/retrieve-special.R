#' Special cases for the magma retrieve function, as laid out in https://mountetna.github.io/magma.html#retrieve
#' @name retrieve_SpecialCases
#' @inheritParams retrieve
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
    projectName,
    token = .get_TOKEN(),
    ...
) {
    
    .retrieve(
        projectName = projectName,
        modelName = "all",
        recordNames = "[]",
        attributeNames = "all",
        format = "json",
        token = token,
        ...)
}

#' @describeIn retrieve_SpecialCases Retrieve the modelNames for a given project
#' @export
retrieveModels <- function(
    projectName,
    token = .get_TOKEN(),
    ...
) {
    names(retrieveTemplate(projectName, token, ...)$models)
}

#' @describeIn retrieve_SpecialCases Retrieve all the identifiers/recordNames for a given project-model pair.
#' @export
retrieveIds <- function(
    projectName,
    modelName,
    token = .get_TOKEN(),
    ...
) {
    
    ### This is very inefficient and should be fixed.
    # This chunk is all that should be required
    .retrieve(
        projectName = projectName,
        modelName = modelName,
        recordNames = "all",
        attributeNames = "identifier",
        format = "tsv",
        token = token,
        ...)[,,drop = TRUE]
}

#' @describeIn retrieve_SpecialCases Retrieve all the attribute options for a given project-model pair.
#' @export
retrieveAttributes <- function(
    projectName,
    modelName,
    token = .get_TOKEN(),
    ...
) {

    rec <- retrieveIds(projectName, modelName, token, ...)[1]
    
    .retrieve(
        projectName, modelName, recordNames = rec,
        attributeNames = "all", format = "tsv",
        names.only = TRUE, token = token, ...)
}
