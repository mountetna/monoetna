#' Assess the downloadable data from magma
#' @name retrieveAttributes
#' @param token your TOKEN from https://janus.ucsf.edu
#' @param projectName the name of the project you would like to download data from, exactly as it appears at https://timur.ucsf.edu
#' @param modelName the name of the piece of the magma model form which to download data
#' @param recordNames,attributeNames particular record(s) or attribute(s) to grab names for
#' @inheritParams retrieve
#' @return a string vector
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token.
#'     retrieveAttributeNames(
#'         projectName = "ipi",
#'         modelName = "patient")
#' }
#' 
#' @importFrom utils URLencode
NULL

#' @describeIn retrieveAttributes Assess the names of attributes downloadable from magma
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

#' @describeIn retrieveAttributes Assess the counts of records for attributes downloadable from magma
#' @export
retrieveAttributeCounts <- function(
    projectName = "ipi",
    modelName = "patient",
    recordNames = "all",
    attributeNames = "all",
    token = .get_TOKEN(),
    ...
) {
    
    if (identical(attributeNames, "all")) {
        attributeNames <- retrieveAttributeNames(
            projectName = projectName, modelName = modelName,
            recordNames = recordNames, token = token, ...)
    }
    
    ### Logic for this function is not complete ###
    lapply(
        attributeNames,
        function(x) {
            dat <- retrieve(
                projectName = projectName, modelName = modelName,
                recordNames = recordNames, attributeNames = x,
                token = token, ...)
            if (is.data.frame(dat)) {
                stop("Function does not work for this data currently")
            } else if (length(attributeNames)==1) {
                sum(unlist(lapply(dat, function(x) length(x)==2)))
            } else {
                stop("Function does not work for this data currently")
            }
        })
}
