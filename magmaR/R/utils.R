.get_TOKEN <- function() {
    # This function asks the user to provide their janus token, and then stores
    # the token as a hidden variable in the workspace.
    # Perhaps a DANGEROUS method, so not necessarily what to do in the end.
    
    if (!exists(".TOKEN")) {
        
        if (interactive()) {
            assign(
                ".TOKEN",
                readline(prompt = "Enter your Janus TOKEN (without quotes):"),
                envir = .GlobalEnv)
        } else {
            stop("Please provide your janus token to the 'token' input or assign it to '.TOKEN'.")
        }
    }
    
    .TOKEN
}

.parse_JSON_tsv <- function(string, names.only = FALSE, connected.only = TRUE) {
    
    rows <- strsplit(string, split = "\n")[[1]]
    
    .split_by_tab <- function(X) {strsplit(X, split = "\t")[[1]]}
    
    cols <- .split_by_tab(rows[1])
    
    if(names.only) {
        return(cols)
    }
    
    ### Parse data
    dat <- sapply(rows[-1], .split_by_tab)
    
    ### Clean data
    if (!is.matrix(dat)) {
        ### When data is not clean, it won't be a matrix.
        dat <- sapply(
            dat,
            function(x) {
                # Add parent attribute = NA for disconnected data
                if (length(x)+1 == length(cols)) {
                    if (connected.only){
                        # Missing entry should be parent. Remove!
                        NULL
                    } else {
                        append(x, values = NA, after = 1)
                    }
                } else {
                    x
                }
            })
        
        if (connected.only) {
            # Remove fully empty records (no parent!)
            dat[sapply(dat, is.null)] <- NULL
            
            # Simplify to a matrix if data is now rectangular.
            dat <- sapply(dat, function(x){x})
        }
    } else {
        if (connected.only) {
            # Remove records with a NULL first column (=parent)
            dat <- dat[!is.null(dat[,1]),]
        }
    }

    ### Add names
    # Hopefully, the data is rectangular
    if (is.matrix(dat)) {
        dat <- as.data.frame(t(dat))
        colnames(dat) <- cols
        rownames(dat) <- seq_len(nrow(dat))
    # "Fallback" options is a list.
    } else if (is.list(dat)) {
        for (i in seq_along(dat)) {
            if (length(dat[[i]]) == length(cols)) {
                names(dat[[i]]) <- cols
            }
        }
        names(dat) <- seq_along(dat)
    }
    
    dat
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
