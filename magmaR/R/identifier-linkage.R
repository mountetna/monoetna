#' @importFrom utils tail

.trace_model_to_proj <- function(
    main_modelName,
    template) {
    
    path <- main_modelName
    while (tail(path,1)!="project") {
        path <- c(path, template$models[[tail(path,1)]]$template$parent)
    }
    
    path
}

.obtain_linkage_paths <- function(
    main_modelName,
    meta_modelName,
    template) {
    
    # Output NA if main_modelName & meta_modelName are the same
    if (main_modelName == meta_modelName) {
        return(NA)
    }
    
    # Get main path
    main_path <- .trace_model_to_proj(main_modelName, template)
    
    # Check for meta model in this path and trim, otherwise, find intersect, output both paths trimmed to intersect.
    if (meta_modelName %in% main_path) {
        ind <- match(meta_modelName, main_path)
        return(list(
            main_path = main_path[seq_len(ind)]
            ))
    } else {
        meta_path <- .trace_model_to_proj(meta_modelName, template)
        
        main_ind <- min(match(meta_path, main_path), na.rm = TRUE)
        meta_ind <- match(main_path[main_ind], meta_path)
        return(list(
            main_path = main_path[seq_len(main_ind)],
            meta_path = meta_path[seq_len(meta_ind)]
            ))
    }
}

#### Next Steps ####

.map_identifiers_by_path <- function(
    path, projectName) {
    
    ids <- query(
            projectName = projectName,
            queryTerms = 
                list(path[1],
                     '::all',
                     path[2],
                     '::identifier'),
            format = "df")
    
    ind <- 2
    
    while (ind < length(path)) {
        
        new_id_map <- query(
            projectName = projectName,
            queryTerms = 
                list(path[ind],
                     '::all',
                     path[ind+1],
                     '::identifier'),
            format = "df")
        
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

# # Obtain data (retrieve or query...)
# 
# .reshape_data_for_multimap_ids <- function(){}

