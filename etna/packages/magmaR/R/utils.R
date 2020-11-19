.get_TOKEN <- function() {
    # This function asks the user to provide their janus token, and then stores
    # the token as a hidden variable in the workspace.
    # Perhaps a DANGEROUS method, so not necessarily what to do in the end.
    
    if (!exists(".MAGMAR_TOKEN")) {
        
        if (interactive()) {
            assign(
                ".MAGMAR_TOKEN",
                readline(prompt = "Enter your Janus TOKEN (without quotes):"),
                envir = .GlobalEnv)
        } else {
            stop("Please provide your janus token to the 'token' input or assign it to '.MAGMAR_TOKEN'.")
        }
    }
    
    .MAGMAR_TOKEN
}

.get_URL <- function() {
    # This function asks the user to provide their janus token, and then stores
    # the token as a hidden variable in the workspace.
    # Perhaps a DANGEROUS method, so not necessarily what to do in the end.
    
    if (!exists(".MAGMAR_URL")) {
        assign(
            ".MAGMAR_URL",
            "https://magma.ucsf.edu",
            envir = .GlobalEnv)
    }
    
    .MAGMAR_URL
}

.perform_curl <- function(
    fxn = c("/retrieve", "/query"),
    requestBody,
    token,
    url.base,
    verbose) {
    
    fxn <- match.arg(fxn)
    
    curl <- RCurl::basicTextGatherer()
    
    RCurl::curlPerform(
        url = paste0(url.base, fxn),
        httpheader = c('Content-Type' = "application/json", 'Authorization' = paste0('Etna ', token)),
        postfields = requestBody,
        writefunction = curl$update,
        verbose = verbose
        )
    
    if (curl$value() == "You are unauthorized") {
        stop("You are unauthorized. If you think this is a mistake, run `rm(.MAGMAR_TOKEN)` or update yout 'token' input, then retry.")
    }
    
    curl
}

.parse_tsv <- function(
    string, names.only = FALSE, connected.only = TRUE) {
    
    ### Parse
    table <- as.data.frame(readr::read_tsv(string))

    if(names.only) {
        return(names(table))
    }
    
    ### Clean data
    # Logic is broken, but tooling is on the way.
    # if (connected.only) {
    #     # Remove if no data for parent.
    #     dat <- dat[!is.null(dat[,1]),]
    # }
    
    table
}

### RENAME
.as_array_unless_all_or_brackets <- function(values) {
    
    if (identical(values, "all")) {
        return("all")
    }
    
    if (identical(values, "[]")) {
        return(I(list()))
    }
    
    if (length(values) > 1) {
        I(values)
    } else {
        I(list(values))
    }
}

### RENAME
.as_array_unless_all_or_identifier <- function(values) {
    
    if (identical(values, "all") || identical(values, "identifier")) {
        return(values)
    }
    
    if (length(values) > 1) {
        I(values)
    } else {
        I(list(values))
    }
}
