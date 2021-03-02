.get_TOKEN <- function() {
    # This function asks the user to provide their janus token, one time, and 
    # then stores the response as a hidden variable in the workspace. For that
    # and future calls, it then provides this hidden variable as its response.
    #
    # Vignettes and error messages make this behavior clear.
    
    if (!exists(".MAGMAR_TOKEN", envir = .GlobalEnv)) {
        
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
    # This function sets up the default magma url as a hidden variable in the
    # workspace. This method of defaulting allows users to adjust their default
    # url by manually setting .MAGMAR_URL <- <different_url> to bypass the
    # behvaior below.
    #
    # Vignettes and error messages make this behavior clear.
    
    if (!exists(".MAGMAR_URL", envir = .GlobalEnv)) {
        assign(
            ".MAGMAR_URL",
            "https://magma.ucsf.edu",
            envir = .GlobalEnv)
    }
    
    .MAGMAR_URL
}

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
    
    if (curl$status_code==401) {
        stop("You are unauthorized. If you think this is a mistake, run `rm(.MAGMAR_TOKEN)` or update your 'token' input, then retry.")
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
