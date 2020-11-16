#' Special cases for the magma retrieve function, as laid out in https://mountetna.github.io/magma.html#retrieve
#' @name retrieve_SpecialCases
#' @inheritParams retrieve
#' @return retrieveTemplate = a list conversion of the project's template json.
#' 
#' retrieveIds = a string vector of record names.
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     
#'     retrieveTemplate(
#'         projectName = "ipi")
#'         
#'     retrieveIds(
#'         projectName = "ipi",
#'         modelName = "rna_seq")
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
        ...)
}

#' @describeIn retrieve_SpecialCases Retrieve all the identifiers for a given project-model pair.
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
        ...)[,,drop = TRUE]
}
