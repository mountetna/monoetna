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


#' (Not complete) Retrieve data from one model of a project that is linked to data from a different, target model of the same project
#' @inheritParams retrieve
#' @param main_modelName,main_recordNames Strings which indicate the target data for which meta-data is desired.
#' @param meta_modelName,meta_attributeNames Strings which indicate the linked "meta"data to retrieve.
#' @param meta_group.by String, an attribute within the meta_model that can be used to group data when there are expected to be more than 1 record of metadata per record of target data.
#' @param ... Additional parameters passed along to the internal `.query()` and `.retrieve()` functions.
#' For troubleshooting or privileged-user purposes only.
#' Options: \code{verbose} (Logical), or \code{url.base} (String which can be user to direct toward production versus staging versus development versions of maagma).
#' @return A dataframe with rows = main_recordNames and columns = either meta_attributeNames or repeats of main_attributeNames to accommodate 1:many mappings.
#' @export
#' @examples
#' 
#' ##### FILLER
#' 
#' # Unless a working TOKEN is hard-coded, or it is in an interactive session,
#' #   this code will not work.
#' 
#' if (interactive()) {
#'     # Running like this will ask for input of your janus token.
#'     retrieveJSON(
#'         projectName = "ipi",
#'         modelName = "patient",
#'         recordNames = "all",
#'         attributeNames = "all",
#'         filter = "")
#' }
#'
retrieveMetadata <- function(
    projectName,
    main_modelName,
    main_recordNames = "all",
    meta_modelName,
    meta_attributeNames = "all",
    meta_group.by = NULL,
    token = .get_TOKEN(),
    ...) {
    
    temp <- retrieveTemplate(projectName, token = token, ...)
    
    # Establish linkage
    paths <- .obtain_linkage_paths(main_modelName, meta_modelName, temp)
    
    meta_id_map <- .map_identifiers_by_path(
        path = paths$meta_path,
        projectName = projectName,
        token = token,
        ...)
    
    main_id_map <- .map_identifiers_by_path(
        path = paths$main_path,
        projectName = projectName,
        token = token,
        ...)
    
    if (!identical(main_recordNames, "all")) {
        main_id_map <- main_id_map[main_id_map[,1] %in% main_recordNames,]
    }
    
    # Subset meta side to matching data
    meta_id_map <- meta_id_map[
        meta_id_map[,ncol(meta_id_map)] %in%
            main_id_map[,ncol(main_id_map)],
    ]
    
    ### Retrieve metadata
    meta_raw <- retrieve(
        projectName = projectName, modelName = meta_modelName,
        recordNames = meta_id_map[,1],
        attributeNames = meta_attributeNames,
        token = token,
        ...)
    
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
    
    ### Expand metadata (columns) per repeated linker ids on the metadata side
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
            # Update column names with _i for all but the main identifier column
            colnames(next_meta) <- paste(colnames(next_meta), i, sep="_")
            colnames(next_meta)[1] <- colnames(meta)[1]
            
            # Add
            meta <- merge(meta, next_meta, all.x = TRUE, by = colnames(meta)[1])
            
            # Update variables for next repeat
            meta_left <- meta_left[-inds,]
            linker_ids_left <- linker_ids_left[-inds]
        }
    }
    
    ### Expand metadata (rows) per repeated ids on the target data side
    output <- merge(main_id_map, meta, all.x = TRUE, by = colnames(meta)[1])
    
    output
}

