#' Special cases for the magma retrieve function as laid out in https://mountetna.github.io/magma.html#retrieve
#' @name retrieve_SpecialCases
#' @param token your TOKEN from https://janus.ucsf.edu
#' @param projectName the name of the project you would like to download data from, exactly as it appears at https://timur.ucsf.edu
#' @param modelName the name of the piece of the magma model form which to download data
#' @inheritParams retrieve
#' @return retrieveTemplate = a jsonlite::fromJSON conversion of template JSON.
#' 
#' retrieveIds = a string vector
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token.
#'     
#'     retrieveTemplate(
#'         projectName = "ipi")
#'         
#'     retrieveIds(
#'         prjectName = "ipi",
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
        ...)
}
