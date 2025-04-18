#' @export
subsetDF_index_targets <- function(data_frame, conditions) {
    # This function converts from the json-list returned by the subsetDataFrame UI component,
    # into the matching indexes of a target dataframe.
    #
    # Inputs:
    # 'conditions' = list structure equivalent of the JSON value output from a vulcan-ui 'SelectionDefinitionPiece'.
    #   If no subsetting was allowed and chosen, it will be an empty list.
    #   Otherwise, it will be an array of targeting constraints where each has:
    #   '$col' = the column name targeted by this constraint.
    #   '$def' = depends on the type of data the targeted column contains, but defines the values to keep
    #        - numeric data: (1) string, "above" for '>' or "exactly" for '>=' for the lower bound; (2) number, lower bound cutoff value; (3) string, "below" for '<' or "exactly" for '<=' for the upper bound; (4) number, lower bound cutoff value.
    #        - string data: (1,2,3,...; however many) string, all elements represent the set of values to be kept.
    #        - boolean data: (1) TRUE or FALSE, whether to keep TRUEs or FALSEs.
    #   '$logic' = NULL for the first contraint, then "AND" or "OR" for all others, which direct how to combine outcomes of the individual constraint definitions.
    #      Per the original design & documentation, the parse ordering is 'all ands first, then all ors', however a potential ToDo (primarily gated behind viable UI-design) would be to give the user more extensive control of logic combination / ordering.
    # 'data_frame' = a data.frame holding the data directly getting subset, or a summary of the ultimate target data which still has structure of..., individual data points in rows and attributes of those data points in columns.
    #
    # Quick Summary:
    #   The function takes in a data_frame describing individual data points, and a set of conditions which define the data points the user wishes to have kept.
    #   It parses the set of data points which meet the user's conditions, and returns the rownames of all data points to retain.  
    #
    # Note: In addition to serving the purpose of subsetting a data.frame directly, the function can also serve as the internal workhorse for subsetting of additional datatypes.
    #   For such cases, a wrapper should be built which:
    #     1. Creates the required 'data_frame' representation where each row equates to a single data point and each column holds the values associated with each data point for feature targetted by 'condition$method's.
    #     2. Call subsetDF_index_targets to identify the data points desired to be kept
    #     3. (Optional, but recommended to be done here rather than leaving to workflow scripts) Perform subsetting of the original data structure based on the "rownames(data_frame)[keep]" return.
    # 
    # Example input for conditions: c(list("col" = "gene1", "def" = c("exactly", 0, "below", 0.02), "logic" = NULL), list("col" = "cell_type", "def" = c("T","B","Myeloid"), "logic" = "AND"))
    # Example2: list() for no subsetting.
    
    types <- c()
    for (i in seq_along(conditions)) {
        types <- c(types, .check_constraint_format(conditions[[i]], i))
    }
    
    # Parse when each condition is True/False.
    raw_calls <- lapply(
        seq_along(conditions),
        function(i) {
            .interpret_constraint(
                conditions[[i]]$col,
                conditions[[i]]$def,
                types[i],
                data_frame
            )
        })

    # Combine all calls by and/or logic
    if (length(conditions)>1) {
        logic <- c("AND", sapply(2:length(conditions), function(i) { conditions[[i]]$logic }))
        ands <- list(1)
        for (i in seq_along(logic)[-1]) {
            if (logic[i]=="AND") {
                ands[[length(ands)]] <- c(ands[[length(ands)]], i)
            } else {
                ands[[length(ands)+1]] <- i
            }
        }
        and_calls <- lapply(ands, function(x) {.map_by_2s(.combine_ands, raw_calls[x])})
        which(.map_by_2s(.combine_ors, and_calls), data_frame)
    } else {
        which(raw_calls[[1]], data_frame)
    }
}

.check_constraint_format <- function(constraint, index) {
    issues <- c()
    type <- "unknown"

    if (!"col" %in% names(constraint)) {
        issues <- c(issues, "'col' missing;")
    }
    if (!"def" %in% names(constraint)) {
        issues <- c(issues, "'def' missing;")
    }
    if (!"logic" %in% names(constraint) || index>1 && !(identical(constraint$logic, "AND") || identical(constraint$logic, "OR"))) {
        issues <- c(issues, "'logic' missing;")
    }

    # End early
    if (length(issues)>0) {
        stop("For index-", index, " constraint:", issues)
    }

    if (any(.is_numeric(constraint$def))) {
        # Numeric data target
        type <- "numeric"
        if (! constraint$def[1] %in% c("exactly", "above")) {
            issues <- c(issues, "numeric 'def' must start with 'exactly' or 'above';")
        }
        if (! constraint$def[3] %in% c("exactly", "below")) {
            issues <- c(issues, "numeric 'def' 3rd element must be 'exactly' or 'below';")
        }
        if (! .is_numeric(constraint$def[2])) {
            issues <- c(issues, "numeric 'def' 2nd element must be numeric;")
        }
        if (! .is_numeric(constraint$def[4])) {
            issues <- c(issues, "numeric 'def' 4th element must be numeric;")
        }
    } else if (all(.is_logical(constraint$def))) {
        # boolean target, actually boolean
        type <- "boolean"
        if (length(constraint$def) != 1) {
            issues <- c(issues, "boolean 'def' must be only one element;")
        }
    } else {
        # string target
        type <- "string"
        if (! all(sapply(constraint$def, is.character)) && length(constraint$def) >= 1) {
            issues <- c(issues, "string 'def' must contain only strings;")
        }
    }

    if (length(issues)>0) {
        stop("For index-", index, " constraint:", issues)
    }

    type
}

.interpret_constraint <- function(col, def, type, data_frame) {

    target_data <- data_frame[,col, drop = TRUE]

    if (all(.is_numeric(target_data))) {
        # Numeric data
        if (type != "numeric") {
            stop("non-numeric constraint 'def' for numeric data,", type)
        }
        target_data <- as.numeric(target_data)

        if (def[1]=="above") {
            left <- target_data > as.numeric(def[2])
        } else {
            left <- target_data >= as.numeric(def[2])
        }

        if (def[3]=="below") {
            right <- target_data < as.numeric(def[4])
        } else {
            right <- target_data <= as.numeric(def[4])
        }

        out <- left & right
    } else if (all(.is_logical(target_data))) {
        # Logical data
        if (type != "boolean") {
            stop("non-boolean constraint 'def' for boolean data")
        }

        match <- FALSE
        if (tolower(def)=="true") {
            match <- TRUE
        }

        if (is.logical(target_data)) {
            out <- target_data==match
        } else {
            out <- tolower(target_data)==tolower(match)
        }
    } else if (type != "string") {
        stop("non-string constraint 'def' for string data")
    } else {
        out <- target_data %in% def
    }

    out & !is.na(out)
}

.combine_ors <- function(clause1, clause2) {
    return(clause1 | clause2)
}
.combine_ands <- function(clause1, clause2) {
    return(clause1 & clause2)
}
.map_by_2s <- function(fxn, calls) {
    out <- calls[[1]]
    i <- 2
    while (i <= length(calls)) {
        out <- fxn(out, calls[[i]])
        i <- i + 1
    }
    return(out)
}

#' @export
subset__scObj <- function(scdata, conditions) {
    # For a SingleCell object, we need to build a data.frame containing the targets of the subsetting conditions.
    # Then we can pass everything along to the normal function.
    # This function makes use of dittoSeq helpers to be viable for scData in Seurat (v2+) and SingleCellExperiment.

    if (any(! c("methods", "logic") %in% names(conditions))){
        stop("Subsetting conditions are not formatted properly.")
    }

    get_data_for_method <- function(method, scdata) {
        target <- method[1]
        if (is.na(target) || is.null(target)) {
            stop("A subsetting condition is incomplete.")
        }

        # Retrieve target data
        assay <- dittoSeq:::.default_assay(scdata)
        if ( !isMeta(target, scdata) && is(scdata, "Seurat") && length(Seurat::Assays(scdata))>1 ) {
            assays <- Seurat::Assays(scdata)
            assay <- assays[which(sapply(assays, function(x) {isGene(target, scdata, assay = x)}))]
        }
        target_data <- dittoSeq:::.var_OR_get_meta_or_gene(target, scdata, assay=assay)
    }

    df <- do.call(data.frame, lapply(
        conditions$methods,
        function(x) {
            get_data_for_method(x, scdata)
        })
    )
    rownames(df) <- colnames(scdata)
    colnames(df) <- vapply(conditions$methods, function(x) { x[1] }, FUN.VALUE = character(1))

    subsetDF_index_targets(df, conditions)
}
