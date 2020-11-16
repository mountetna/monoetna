#' Determine the attribute options associated with a given project-model pair.
#' @inheritParams retrieve
#' @return a string vector
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     retrieveAttributeNames(
#'         projectName = "ipi",
#'         modelName = "patient")
#' }
#' 
#' @export
retrieveAttributeNames <- function(
    projectName = "ipi",
    modelName = "patient",
    recordNames = "all",
    token = .get_TOKEN(),
    ...
) {

    .retrieve(
        projectName, modelName, recordNames,
        attributeNames = "all", format = "tsv",
        names.only = TRUE, token = token, ...)
}
