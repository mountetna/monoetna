#' Retrieve data from one model of a project, metadata, for records linked to data from a different, target model of the same project
#' @inheritParams retrieve
#' @param target_modelName,target_recordNames Strings which indicate the 'target' data that meta-data is desired to be linked to.
#' @param meta_modelName,meta_attributeNames Strings which indicate the linked "meta"data to retrieve.
#' @param meta_group.by \emph{Not implemented yet.} String, an attribute within the meta_model that can be used to group data when there are expected to be more than 1 record of metadata per record of target data.
#' @param ... Additional parameters passed along to the internal `.query()` and `.retrieve()` functions.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{verbose} (Logical), or \code{url.base} (String which can be user to direct toward production versus staging versus development versions of magma).
#' @return A dataframe with rows = target_recordNames and columns = either meta_attributeNames or repeats of meta_attributeNames to accommodate 1:many mappings of target data records to metadata records.
#' @export
#' @importFrom utils tail
#' @examples
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token one time.
#'     retrieveMetadata(
#'         projectName = "ipi",
#'         meta_modelName = "sample",
#'         meta_attributeNames = "tissue_type",
#'         target_modelName = "rna_seq",
#'         filter = "")
#' }
#'
retrieveMetadata <- function(
    projectName,
    meta_modelName,
    meta_attributeNames = "all",
    target_modelName,
    target_recordNames = "all",
    meta_group.by = NULL,
    token = .get_TOKEN(),
    ...) {
    
    temp <- retrieveTemplate(projectName, token = token, ...)
    
    ### Establish linkage and retrieve identifier mappings
    paths <- .obtain_linkage_paths(target_modelName, meta_modelName, temp)
    separate_branches <- length(paths) == 2
    
    ### Target data
    target_id_map <- .map_identifiers_by_path(
        path = paths$target_path,
        projectName = projectName,
        token = token,
        ...)
    
    if (!identical(target_recordNames, "all")) {
        target_id_map <- target_id_map[target_id_map[,1] %in% target_recordNames,]
    }
    
    ### Metadata
    # Determine linked records
    if (separate_branches) {
        meta_id_map <- .map_identifiers_by_path(
        path = paths$meta_path,
        projectName = projectName,
        token = token,
        ...)
    
        # Subset meta side to matching data
        meta_id_map <- meta_id_map[
            meta_id_map[,ncol(meta_id_map)] %in%
                target_id_map[,ncol(target_id_map)],
        ]
        
        meta_record_names <- meta_id_map[,1, drop = TRUE]
        
        # Create combined id_map with the meta columns at end, in reverse order.
        target_id_map <- merge(
            target_id_map, meta_id_map[,rev(paths$meta_path)],
            by = colnames(target_id_map)[ncol(target_id_map)])
    } else {
        
        # The target meta_model is at the top of the target_id_map.
        meta_record_names <- target_id_map[,ncol(target_id_map), drop = TRUE]
        meta_record_names <- meta_record_names[!is.na(meta_record_names)]
    }
    
    ### Retrieve metadata
    meta_raw <- retrieve(
        projectName = projectName, modelName = meta_modelName,
        recordNames = meta_record_names,
        attributeNames = meta_attributeNames,
        token = token,
        ...)
    # Find identifier column
    meta_id_col_name <- temp$models[[meta_modelName]]$template$identifier
    id_col_index <- match(meta_id_col_name, colnames(meta_raw))
    
    ### Ensure metadata has 1 row per linked target data row
    meta <- if (separate_branches) {
        .expand_metadata_to_have_1row_per_id(meta_raw, meta_id_map, id_col_index)
    } else {
        meta_raw
    }
    
    ### Output in order that matches target_recordNames order
    
    # Ensure no duplicated column names, except for the matching id column.
    meta <- .trim_duplicated_columns(meta, target_id_map, meta_id_col_name)
    
    # merge with order coming from the the target_id_map
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
    path, projectName, token = .get_TOKEN(), ...) {
    
    ids <- query(
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

.expand_metadata_to_have_1row_per_id <- function(meta_raw, meta_id_map, id_col_index){
    ### Expand metadata (columns direction) until 1 row per id.
    
    ## Method: By a user-provided grouping 
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
    
    ## Method: Simply expand rightward by shifting all data for repeated ids to new columns.
    linker_ids <- meta_id_map[,ncol(meta_id_map)]
    if (any(duplicated(linker_ids))) {
        
        inds <- match(unique(linker_ids), linker_ids)
        
        meta <- meta_raw[inds,]
        
        meta_left <- meta_raw[-inds,]
        linker_ids_left <- linker_ids[-inds]
        
        for (i in 2:max(table(linker_ids))) {
            
            inds <- match(unique(linker_ids_left), linker_ids_left)
            
            # Extract next chunk
            next_meta <- meta_left[inds,]
            # Update column names with _i for all but the column containing Ids
            colnames(next_meta) <- paste(colnames(next_meta), i, sep="_")
            colnames(next_meta)[id_col_index] <- colnames(meta)[id_col_index]
            
            # Add
            meta <- merge(meta, next_meta, all.x = TRUE, by = colnames(meta)[id_col_index])
            
            # Update variables for next repeat
            meta_left <- meta_left[-inds,]
            linker_ids_left <- linker_ids_left[-inds]
        }
    } else {
        meta <- meta_raw
    }
    
    meta
}

.trim_duplicated_columns <- function(meta, target_id_map, meta_id_col_name) {
    
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
