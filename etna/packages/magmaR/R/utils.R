.perform_curl_get <- function(
    fxn = c("/retrieve", "/query", "/update", "/projects"),
    target,
    requestBody,
    parse = TRUE,
    verbose = FALSE) {
    
    fxn <- match.arg(fxn)
    
    opts <- target$opts
    opts$postfields <- requestBody
    
    # Set
    curl <- crul::HttpClient$new(
        url = target$url,
        headers = list(
            'Content-Type' = "application/json",
            'Authorization' = paste0('Etna ', target$token)),
        opts = opts
        )
    
    # Perform
    curl <- curl$get(path = fxn)
    
    # Summarize
    if (verbose) {
        if (curl$success()) {
            cat(paste0(fxn, ": successful."))
        } else {
            cat(fxn, ":\n")
            print(curl$status_http())
        }
    }
    
    if (curl$status_code %in% c(302,401)) {
        stop("You are unauthorized. Update your 'token' input with 'magmaRset()', then retry.")
    }
    
    # Parse
    if (parse) {
        output <- curl$parse(encoding = "UTF-8")
    } else {
        output <- curl
    }
    
    output
}

.parse_tsv <- function(
    string, names.only = FALSE) {
    
    # Parse
    table <- read.csv(text = string, sep = "\t")

    if(names.only) {
        return(names(table))
    }
    
    table
}

.get_sysenv_or_mock <- function(target) {
    ifelse(
        identical(Sys.getenv(target),""),
        paste0("fake", tolower(target)),
        Sys.getenv(target))
}
