#' Analogous to the '/update' function of magma
#' @description Analogous to the '/update' function of magma, allows data to be sent to magma (by users with at least "editor" authorization).
#' @inheritParams updateMatrix
#' @inheritParams retrieve
#' @param revisions A list of named lists containing the data to be updated.
#' 
#' List structure:
#' \itemize{
#' \item top level name(s): modelNames, can be 1 or more.
#' \item 2nd level name(s): recordNames, can be 1 or more.
#' \item 3rd level name(s) & contents: the attributes to update & the values to use.
#' }
#' 
#' See \url{https://mountetna.github.io/magma.html#update} for additional formatting details.
#' 
#' @param auto.proceed Logical. When set to TRUE, the function does not ask before proceeding forward with the 'magma/update'.
#' @return None directly.
#' 
#' The function sends data to magma, and the only outputs are information reported via the console.
#' @details This function mimics the activity of the magma/update function, documented here \url{https://mountetna.github.io/magma.html#update},
#' with the main difference being that the \code{revisions} input should be in nested list format rather than nested hash (because R does not support hash structures).
#' 
#' Internally, the function:
#' 
#' 1. Summarizes records of each model that will be targeted for updating.
#' 
#' 2. Prompts the user before proceeding (unless \code{auto.proceed} is set to \code{TRUE})
#' 
#' 3. Directly passes its inputs along to magma/update via a curl request.
#' 
#' @seealso
#' \url{https://mountetna.github.io/magma.html#update} for documentation of the underlying magma/update function.
#' 
#' \code{\link{updateMatrix}} for a matrix-dedicated version of this function which can be provided a matrix, or matrix's file location, directly.
#' @export
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     # Note that you likely do not have write-permissions for the 'example'
#'     # project, so this code can be expected to give an authorization error.
#' 
#'     updateValues(
#'         target = magma,
#'         projectName = "example",
#'         revisions = list(
#'             # model
#'             'rna_seq' = list(
#'                 # record
#'                 'EXAMPLE-HS1-WB1-RSQ1' = list(
#'                     # attribute
#'                     'fraction' = list(
#'                         # value(s)
#'                         "Tcells"
#'                         )
#'                     )
#'                 )
#'             )
#'         )
#' }
#'
updateValues <- function(
    target,
    projectName,
    revisions = list(),
    auto.proceed = FALSE,
    ...
) {
    
    ### Summarize per model
    lapply(names(revisions), function(model) {
        
        current_ids <- retrieveIds(
            target, projectName, model)
        
        .summarize_model_values(revisions[[model]], model, current_ids)
    })
    
    ### Check if should move forward
    go <- .ask_before_proceeding(auto.proceed)
    if (! tolower(go) %in% c("y", "yes")) {
        return("No /update performed.")
    }
    
    .update(
        target = target,
        projectName = projectName,
        revisions = revisions,
        ...)
}

.summarize_model_values <- function(model_revs, modelName, model_ids) {
    
    ### Report how many records will be updated
    num_recs <- num_current_recs <- length(model_revs)
    rec_names <- rec_names_current <- names(model_revs)
    
    # Check if any would be be new
    if (! all(rec_names %in% model_ids)) {
        
        rec_names_current <- rec_names[rec_names %in% model_ids]
        
        rec_names_new <- rec_names[!rec_names %in% model_ids]
        num_new <- length(rec_names_new)
        
        ### Summarize for NEW records
        cat("For model \"", modelName, "\", this update() will create (or update) ", num_new, " NEW (or orphan) records:\n    ",
            paste0(rec_names_new, collapse = "\n    "),
            "\nWARNING: Check the above carefully. Once created, there is no easy way to remove records from magma.\n",
            sep="")
        
        num_current_recs <- num_recs - num_new
    }
    
    ### Summarize for current records.
    cat("For model \"", modelName, "\", this update() will update ", num_current_recs, " records",
        if (num_current_recs==0) {
                ".\n"
            } else {
                paste0(":\n    ", paste0(rec_names_current, collapse = "\n    "))
            },
        "\n",
        sep="")
    
    model_revs
}

.ask_before_proceeding <- function(auto.proceed = FALSE) {
    
    if (!auto.proceed) {
        
        if (!interactive()) {
            warning("To run in non-interactive mode, set 'auto.proceed = TRUE'.")
        }
        
        return(readline("Proceed, Y/n? "))
    } else {
        return("Y")
    }
}

.update <- function(
    target,
    projectName,
    revisions,
    request.only = FALSE,
    json.params.only = FALSE,
    verbose = TRUE
) {
    
    ### Put together the request 
    jsonParams <- list(
        project_name = projectName,
        revisions = revisions)
    
    requestBody <- jsonlite::toJSON(jsonParams, auto_unbox = TRUE)
    
    ### Output here if requested.
    if (request.only) {
        return(requestBody)
    }
    if (json.params.only) {
        return(jsonParams)
    }
    
    ### Perform '\update'
    curl_out <- .perform_curl_get(
        fxn = "/update",
        target, requestBody,
        parse = FALSE, verbose = verbose)
}
