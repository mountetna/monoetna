#' Helper function that retrieves all the projectName options which a user has access to, from janus.
#' @inheritParams retrieve
#' @param all_cols Logical. Sets whether to report back the entire table instead of the default behavior where the column 'cc_text' is removed for brevity.
#' @param verbose Logical. Sets whether to report the status of the '/api/user/projects' curl request sent to janus.
#' @details
#' This function takes in the user's \code{target} containing their authorization token, and a url targeting either magma or janus.
#' It then converts the given url to target janus, and makes a curl request to <janus-url>/api/user/projects in order to return which projects a user can access.
#' @return A data.frame where elements of the 'project_name' column reflect what can be given to \code{projectName} inputs of other magmaR functions.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     retrieveProjects(target = magmaRset())
#' }
#' 
retrieveProjects <- function(
    target,
    all_cols = FALSE,
    verbose = FALSE
) {
    
    # Swap in 'janus' for 'magma' in the target url
    target$url <- gsub("magma", "janus", target$url)
    
    return <- .perform_curl_get(
        fxn = c("/api/user/projects"),
        target = target,
        requestBody = NULL,
        parse = TRUE,
        verbose = verbose)
    
    df <- jsonlite::fromJSON(
        return
        )$projects

    if (all_cols) {
        df
    } else {
        df[,!colnames(df) %in% c('cc_text')]
    }
}
