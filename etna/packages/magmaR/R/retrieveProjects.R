#' Retrieve all the projectName options that a user has access to, from janus.
#' @inheritParams retrieve
#' @param url.base Single string. A url which directs towards production versus staging versus development etna environments. See \code{\link{magma-environments}}.
#' 
#' Note that here, we actually pull from 'janus' rather than 'magma', but the function will correct any magma-targeting urls to target janus.
#' @param verbose Logical. Sets whether to report the status of the '/projects' curl request sent to janus.
#' @details
#' This function takes in a user's authorization token, and a url targetting either magma or janus.
#' It then converts that url to target janus, and makes a curl request to <janus-url>/projects in order to return which projects a user can access.
#' @return A data.frame where elements of the 'project_name' column reflect what can be given to \code{projectName} inputs of other magmaR functions.
#' @export
#' @examples
#' 
#' if (interactive()) {
#' 
#'     # Running like this will ask for input of your janus token one time.
#'     retrieveProjects()
#' }
#' 
retrieveProjects <- function(
    token = .get_TOKEN(),
    url.base = .get_URL(),
    verbose = FALSE
) {
    
    # Swap in 'janus' for 'magma' in the target url
    url.base <- gsub("magma", "janus", url.base)
    
    return <- .perform_curl_get(
        fxn = c("/projects"),
        requestBody = NULL,
        token = token,
        url.base = url.base,
        parse = TRUE,
        verbose = verbose)
    
    jsonlite::fromJSON(
        return
        )$projects
}
