.jsonify <- function(list) {
    jsonlite::toJSON(list, auto_unbox = TRUE, null = 'null', na = 'string')
}

.match_expected_recName_structure <- function(values) {

    if (identical(values,"[]")) {
        return(I(list()))
    }

    .match_expected_common(values)
}

.match_expected_attName_structure <- function(values) {
    # "identifier" -> stays this way
    # anything else -> ["value1", "value2", etc.]

    if (identical(values, "identifier")) {
        return(values)
    }

    .match_expected_common(values)
}

.match_expected_common <- function(values) {
    # all -> all
    # vector -> I(vector)
    # single -> I(list(single))

    if (identical(values, "all")) {
        return(values)
    }

    if (length(values) > 1) {
        I(values)
    } else {
        I(list(values))
    }
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

.is_table_model <- function(
    target,
    projectName,
    modelName,
    template = NULL
) {

    if (modelName=="project") {
        return(FALSE)
    }

    if (identical(template, NULL)) {
        template <- retrieveTemplate(target, projectName)
    }

    if (!modelName %in% names(template$models)) {
        stop("'", modelName, "' is not a model of the '", projectName, "' project.")
    }
    parentModelName <- template$models[[modelName]]$template$parent
    template$models[[parentModelName]]$template$attributes[[modelName]]$attribute_type=='table'
}

# Transform into the nested list format
# Note: Do not supply recordNames directly to vapply as any "-" will be
#   converted to "."
.df_to_revisions <- function(DF, modelName) {

    DF_noID <- DF[, seq_len(ncol(DF))[-1], drop = FALSE]
    # For each row of the DataFrame...
    recs <- lapply(
        seq_len(nrow(DF_noID)),
        function(x) {
            # Make the contents of cols 2:end a list of attribute values, and for each
            # attribute value slot, make it a list if length is >1.
            atts <- lapply(
                seq_len(ncol(DF_noID)),
                function(y) {
                    DF_noID[x,y]
                })
            names(atts) <- colnames(DF_noID)
            atts
        })
    names(recs) <- DF[,1, drop = TRUE]

    setNames(list(recs), modelName)
}

.get_sysenv_or_mock <- function(target) {
    ifelse(
        identical(Sys.getenv(target),""),
        paste0("fake", tolower(target)),
        Sys.getenv(target))
}
