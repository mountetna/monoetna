anyAPI <- function(
        app_url,
        api_path,
        target,
        method = c('get', 'post', 'delete'),
        requestBody = NULL,
        raw = FALSE,
        jsonparse = TRUE,
        verbose = FALSE) {

    method <- match.arg(method)

    opts <- target$opts
    opts$postfields <- requestBody

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
    curl <- curl$get(path = api_path)

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
