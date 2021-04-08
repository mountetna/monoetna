#' Download data from magma of one model, but transformed into the shape of a different model's records.
#' @description Retrieve data from one model ("meta") transformed into the shape of linked records of a different model ("target").
#' For example, one could get subject-level information for an RNAseq counts matrix with this function.
#' The output would contain columns of subject-level attributes, and rows that are the RNAseq-model records.
#' @inheritParams retrieve
#' @param target_modelName,target_recordNames Strings which indicate the "target" data that meta-data is desired to be reshaped into.
#' They work the same as inputs of other functions without the \code{target_} portion, and
#' these inputs ultimately set which records of "meta" model data to actually obtain.
#' @param meta_modelName,meta_attributeNames Strings which indicate the "meta"data to retrieve.
#' They work the same as inputs of other functions without the \code{meta_} portion.
#' @details This function retrieves data from one model (the "meta" model) transformed so that rows of the returned dataframe correspond to records of a different model (the "target" model).
#' 
#' Internally, it first determines the path, through child-parent model linkages, to navigate from the meta_model to the target_model.
#' 
#' Then, it performs calls to \code{\link{query}} in order to retrieve identifier linkage along that path.
#' The identifier linkage is turned into a dataframe of identifier traces.
#' 
#' Next, it performs a call to \code{\link{retrieve}} to obtain the wanted metadata as a dataframe.
#' 
#' (If linkage paths would create any 1:many mappings of target data records to metadata records, data of extra records are shifted "rightwards" into columns appended with "_#" in their names.
#' This is a reliable, though imperfect, method so we hope to implement alternatives in the future.)
#' 
#' Finally, the dataframe of linkage path identifiers is merged with the metadata dataframe, reshaping the metadata to properly have one row per requested \code{target_recordName}.
#' @return A dataframe with rows = \code{target_recordNames} and columns = model identifiers and either \code{meta_attributeNames} or repeats of \code{meta_attributeNames_#} when there are 1:many mappings of target data records to metadata records.
#' @seealso
#' \code{\link{retrieve}} and \code{\link{retrieveMatrix}} which will likely be useful for retrieving associated "target" data.
#' 
#' @export
#' @importFrom utils tail
#' @examples
#' 
#' if (interactive()) {
#'     # First, we use magmaRset to create an object which will tell other magmaR
#'     #  functions our authentication token (as well as some other optional bits).
#'     # When run in this way, it will ask you to give your token.
#'     magma <- magmaRset()
#'     
#'     # Running like this will ask for input of your janus token one time.
#'     retrieveMetadata(
#'         target = magma,
#'         projectName = "example",
#'         meta_modelName = "subject",
#'         meta_attributeNames = "group",
#'         target_modelName = "rna_seq",
#'         target_recordNames = "all")
#' }
#'
retrieveMetadata <- function(
    target,
    projectName,
    meta_modelName,
    meta_attributeNames = "all",
    target_modelName,
    target_recordNames = "all",
    ...) {
    
    temp <- retrieveTemplate(target, projectName)
    
    ### Establish linkage and retrieve identifier mappings
    paths <- .obtain_linkage_paths(target_modelName, meta_modelName, temp)
    separate_branches <- length(paths) == 2
    
    ### Target data
    target_id_map <- .map_identifiers_by_path(
        target = target,
        path = paths$target_path,
        projectName = projectName,
        ...)
    
    if (!identical(target_recordNames, "all")) {
        target_id_map <- target_id_map[target_id_map[,1] %in% target_recordNames,]
    }
    
    ### Metadata
    # Determine linked records
    if (separate_branches) {
        meta_id_map <- .map_identifiers_by_path(
            target = target,
            path = paths$meta_path,
            projectName = projectName,
            ...)
    
        # Subset meta side to matching data
        meta_id_map <- meta_id_map[
            meta_id_map[,ncol(meta_id_map)] %in%
                target_id_map[,ncol(target_id_map)],
        ]
        
        meta_record_names <- meta_id_map[,1, drop = TRUE]
    } else {
        
        # The target meta_model is at the top of the target_id_map.
        meta_record_names <- target_id_map[,ncol(target_id_map), drop = TRUE]
        meta_record_names <- meta_record_names[!is.na(meta_record_names)]
    }
    
    ### Retrieve metadata
    meta_raw <- retrieve(
        target = target, projectName = projectName, modelName = meta_modelName,
        recordNames = meta_record_names,
        attributeNames = meta_attributeNames)
    
    ### Ensure metadata has 1 row per linked target data row
    # Find identifier column
    meta_id_col_name <- temp$models[[meta_modelName]]$template$identifier
    id_col_index <- match(meta_id_col_name, colnames(meta_raw))
    # Expand
    meta <- if (separate_branches) {
        .expand_metadata_to_have_1row_per_id(meta_raw, meta_id_map, id_col_index)
    } else {
        meta_raw
    }
    
    ### Output in order that matches target_recordNames order
    
    if (separate_branches) {
        # merge meta df with meta_id_map
        # after ensuring no duplicated column names, except for the id column.
        meta <- .trim_duplicated_columns(meta, meta_id_map, meta_id_col_name)
        meta <-merge(
            meta_id_map, meta,
            by.x = colnames(meta_id_map)[1],
            all.y = TRUE,
            by.y = meta_id_col_name)
        meta_id_col_name <- colnames(meta_id_map)[ncol(meta_id_map)]
    }
    
    # merge with order coming from the the target_id_map
    # after ensuring no duplicated column names, except for the id column.
    meta <- .trim_duplicated_columns(meta, target_id_map, meta_id_col_name)
    output <- merge(
        target_id_map, meta,
        all.x = TRUE,
        by.x = colnames(target_id_map)[ncol(target_id_map)],
        by.y = meta_id_col_name)
    
    output
}

.trace_model_to_proj <- function(
    target_modelName,
    template) {
    
    path <- target_modelName
    while (tail(path,1)!="project") {
        path <- c(path, template$models[[tail(path,1)]]$template$parent)
    }
    
    path
}

.obtain_linkage_paths <- function(
    target_modelName,
    meta_modelName,
    template) {
    
    # Output NA if target_modelName & meta_modelName are the same
    if (target_modelName == meta_modelName) {
        return(NA)
    }
    
    # Get target path
    target_path <- .trace_model_to_proj(target_modelName, template)
    
    # Check for meta model in this path and trim, otherwise, find intersect, output both paths trimmed to intersect.
    if (meta_modelName %in% target_path) {
        ind <- match(meta_modelName, target_path)
        return(list(
            target_path = target_path[seq_len(ind)]
            ))
    } else {
        meta_path <- .trace_model_to_proj(meta_modelName, template)
        
        target_ind <- min(match(meta_path, target_path), na.rm = TRUE)
        meta_ind <- match(target_path[target_ind], meta_path)
        return(list(
            target_path = target_path[seq_len(target_ind)],
            meta_path = meta_path[seq_len(meta_ind)]
            ))
    }
}

.map_identifiers_by_path <- function(
    target, path, projectName, ...) {
    # Function takes in a vector of model names, 'path', for a given project,
    # then makes query() calls to obtain child-parent identifier linkage for
    # all successive pairs of target models. Output is a data.frame which can
    # serve as a dictionary of the linkage between recordNames of any of the
    # models named in 'paths'
    
    ids <- query(
        target = target,
        projectName = projectName,
        queryTerms = 
            list(path[1],
                 '::all',
                 path[2],
                 '::identifier'),
        format = "df",
        ...)
    
    ind <- 2
    
    while (ind < length(path)) {
        
        new_id_map <- query(
            target = target,
            projectName = projectName,
            queryTerms = 
                list(path[ind],
                     '::all',
                     path[ind+1],
                     '::identifier'),
            format = "df",
            ...)
        
        if (any(duplicated(new_id_map[1,]))) {
            stop("Algorithm issue: mappings going upwards not unique")
        }
        
        # Reorder by previous, and add NAs for missing data
        order <- match(ids[,2], new_id_map[,1])
        order <- order[!is.na(order)]
        
        new_ids <- array(NA, dim = nrow(ids))
        new_ids[ids[,2] %in% new_id_map[,1]] <- new_id_map[order ,2]
        
        
        ids <- cbind(ids, new_ids)
        ind <- ind + 1
    }
    
    colnames(ids) <- path
    
    ids
}

.expand_metadata_to_have_1row_per_id <- function(meta_raw, meta_id_map, id_col_index) {
    ### Expand metadata (columns direction) until 1 row per id.
    
    ## Method: By a user-provided grouping 
    # NOT IMPLEMENTED
    # Retrieve grouping data if not already in 'meta'
    # grouping <- if (meta_group.by %in% colnames(meta)) {
    #     meta[,meta_group.by]
    # } else {
    #     retrieve(
    #         projectName = projectName, modelName = meta_modelName,
    #         recordNames = meta_id_map[,1],
    #         attributeNames = meta_attributeNames,
    #         ...
    #         )[,meta_group.by]
    # }
    
    ## Method: Simply expand rightward by shifting all data for repeated ids to
    # new columns.
    
    linker_ids <- meta_id_map[,ncol(meta_id_map)]
    
    if (any(duplicated(linker_ids))) {
        
        # Initialize: Determine first current rows per target final rows, and collect
        inds <- match(unique(linker_ids), linker_ids)
        meta <- meta_raw[inds,]
        
        # Update knowledge of what remains to be added
        meta_left <- meta_raw[-inds,]
        linker_ids_left <- linker_ids[-inds]
        
        for (i in 2:max(table(linker_ids))) {
            
            # Determine next current rows per target final rows
            inds <- match(unique(linker_ids_left), linker_ids_left)
            
            # Collect next chunk
            next_meta <- meta_left[inds,]
            
            # Update column names with _i for all but the column containing Ids
            colnames(next_meta) <- paste(colnames(next_meta), i, sep="_")
            colnames(next_meta)[id_col_index] <- colnames(meta)[id_col_index]
            
            # Add
            meta <- merge(meta, next_meta, all.x = TRUE, by = colnames(meta)[id_col_index])
            
            # Update knowledge of what remains to be added
            meta_left <- meta_left[-inds,]
            linker_ids_left <- linker_ids_left[-inds]
        }
    } else {
        
        meta <- meta_raw
    }
    
    meta
}

.trim_duplicated_columns <- function(meta, target_id_map, meta_id_col_name) {
    # Takes two data.frames, 'meta' & 'target_id_map', and a column name to
    # ignore, 'meta_id_col_name'. It trims any columns from meta, other than
    # the ignored column, that have the same name as a column of target_id_map.
    
    potential_columns <- c(colnames(meta), colnames(target_id_map))
    
    if (any(duplicated(potential_columns))) {
        
        ind <- which(duplicated(potential_columns, fromLast = TRUE))
        # Remove the id column index, from this set.
        ind <- setdiff(ind, grep(meta_id_col_name, colnames(meta)))
        
        # Exclude the ind columns from the meta df, if needed
        if (length(ind)>0) {
            meta <- meta[,-ind]
        }
    }
    
    meta
}
