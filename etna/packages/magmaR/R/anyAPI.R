#' (Advanced) A generalized function for interacting directly with any Data Library App API.
#' @param appName String providing the name of the app with the API you wish to target
#' @param api String providing the endpoint to hit.  e.g. "retrieve" or "gnomon/<project>/rules"
#' @param method "get", "post", or "delete". String providing the http protocol used for this API
#' @param requestBody JSON-format String OR a named list providing any additional parameters used by the given API.
#' When given as a list, where names are the parameter names and contents are the intended values, the contents will be converted to a JSON string internally.
#' @param jsonparse Logical setting whether the return should be converted from JSON-string format into an analogous R format
#' @inheritParams retrieve
#' @return Default = Any structure determined by the behavior of the given API, OR a crul response object if \code{raw} is set to \code{TRUE}.
#' @examples
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'
#'     # Retrieving gnomon rules of a project via the gnomon#rules API
#'     rules_raw <- anyAPI("magma", paste0("gnomon/example/rules"), prod, 'get')
#'     rules <- rules_raw$rules
#'     print(rules)
#' }
#' @export
anyAPI <- function(
        appName = c('magma', 'janus', 'metis', 'timur', 'vulcan', 'polyphemus', 'vesta', 'etna'),
        api,
        target,
        method = c('get', 'post', 'delete'),
        requestBody = NULL,
        raw = FALSE,
        jsonparse = TRUE,
        verbose = FALSE) {

    appName <- match.arg(appName)
    method <- match.arg(method)

    app_url <- gsub('magma', appName, target$url)

    opts <- target$opts

    opts$postfields <- if (is.list(requestBody)) {
        .jsonify(requestBody)
    } else {
        requestBody
    }

    # Set
    curl <- crul::HttpClient$new(
        url = app_url,
        headers = list(
            'Content-Type' = "application/json",
            'Authorization' = paste0('Etna ', target$token)),
        opts = opts
    )

    # Perform
    fxn <- switch(
        method,
        get = curl$get,
        post = curl$post,
        delete = curl$delete)
    curl <- fxn(path = api)

    # Summarize
    if (verbose) {
        if (curl$success()) {
            cat(paste0(url, ": successful."))
        } else {
            cat(url, ":\n")
            print(curl$status_http())
        }
    }

    if (curl$status_code %in% c(302,401)) {
        stop("You are unauthorized. Update your 'token' input with 'magmaRset()', then retry.")
    }

    # Parse
    if (!raw) {
        output <- curl$parse(encoding = "UTF-8")
        if (jsonparse) {
            output <- jsonlite::fromJSON(output)
        }
    } else {
        output <- curl
    }

    output
}
