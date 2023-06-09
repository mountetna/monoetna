subsetDF_index_targets <- function(data_frame, conditions) {
    # This function converts from the json-list returned by the subsetDataFrame UI component,
    # into the matching indexes of a target dataframe.
    #
    # Inputs:
    # 'conditions' = list structure equivalent of the JSON value output from a vulcan-ui 'subsetDataFramePiece'.
    #   When no subsetting was wanted by the user, it will be an empty list. 
    #   Otherwise, it will contain 'methods' and sometimes 'logic'.
    #   'conditions$methods' = list of individual subsetting contraint definitions.
    #      All definitions start with (1) a string representing the name of a column 'data_frame' that they target.
    #      Remaining elements of a definition depend on the type of data the targetted column contains:
    #        - numeric data: (2) string, "above" for '>' or "exactly" for '>=' for the lower bound; (3) number, lower bound cutoff value; (4) string, "below" for '<' or "exactly" for '<=' for the upper bound; (5) number, lower bound cutoff value.
    #        - string data: (2,3,4,...; however many) string, all elements represent the set of values to be kept.
    #        - boolean data: (2) TRUE or FALSE, whether to keep TRUEs or FALSEs.
    #   'conditions$logic' = string vector of "and"s and/or "or"s, which directs how to combine outcomes of individual 'method' definitions and which will only be present if there are at least 2 'methods'.
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
    #     3. (Optional, but recommended to be done here rather than leaving to workflow scripts) Perform subsetting of the original datat structure based on the "rownames(data_frame)[keep]" return.
    # 
    # Example input for conditions: list(methods = list(c("groups", "g1"), c("ACAP1", "exactly", "0", "below", "1"), c("RNA_snn_res.1", "0", "1")), logic = c("and", "and"))
    # Example2: list() for no subsetting.
    
    if (any(! c("methods", "logic") %in% names(conditions))){
        stop("Subsetting conditions are not formatted properly.")
    }
    
    # Parse when each condition is True/False.
    interpret_method <- function(method) {
        if (is.na(method[1]) || is.null(method[1])) stop("A subsetting condition is incomplete.")
        target_data <- data_frame[,method[1], drop = TRUE]
        # Numeric data
        if (is.numeric(target_data)) {
            if ((length(method)!=5)) {
                stop("A subsetting condition for numeric-type data is not formatted properly.")
            }

            if (method[2]=="above") {
                left <- target_data > as.numeric(method[3])
            } else {
                left <- target_data >= as.numeric(method[3])
            }

            if (method[4]=="below") {
                right <- target_data < as.numeric(method[5])
            } else {
                right <- target_data <= as.numeric(method[5])
            }

            return(left & right)
        }
        # Logical data
        if (is.logical(target_data)) {
            if (length(method)!=2 || !tolower(method[2]) %in% c("true", "false")) {
                stop("A subsetting condition for logical-type data was left incomplete.")
            }

            match <- FALSE
            if (tolower(method[2])=="true") {
                match <- TRUE
            }

            return(target_data==match)
        }
        # Remainder = String data
        if (length(method)<2) {
            stop("A subsetting condition for string-type data was left incomplete, or a condition targetting an unimplemented data type was left in place.")
        }
        matches <- method[-1]
        return(target_data %in% matches)
    }
    raw_calls <- lapply(conditions$methods, interpret_method)
    
    # Combine all calls by and/or logic
    combine_ors <- function(clause1, clause2) {
        return(clause1 | clause2)
    }
    combine_ands <- function(clause1, clause2) {
        return(clause1 & clause2)
    }   
    map_by_2s <- function(fxn, calls) {
        out <- calls[[1]]
        i <- 2
        while (i <= length(calls)) {
            out <- fxn(out, calls[[i]])
            i <- i + 1
        }
        return(out)
    }
             
    logic <- conditions$logic
    ands <- list(1)
    if (!is.list(logic)) { # Formatting: list(list())
        for (i in seq_along(logic)) {
            if (logic[i]=="and") {
                ands[[length(ands)]] <- c(ands[[length(ands)]], i)
            } else {
                ands[[length(ands)+1]] <- i
            }
        }
    }
    and_calls <- lapply(ands, function(x) {map_by_2s(combine_ands, raw_calls[x])})

    dittoViz:::.which_rows(map_by_2s(combine_ors, and_calls), data_frame)
}

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
