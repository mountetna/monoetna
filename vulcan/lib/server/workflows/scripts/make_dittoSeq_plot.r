library(Seurat)
library(dataflow)
library(dittoSeq)
library(plotly)

########## Should go into archimedes-r
subsetDF_index_targets <- function(data_frame, conditions) {
    # '''
    # This function converts from the json-list returned by the subsetDataFrame UI component,
    # into the matching indexes of the target dataframe.
    # '''
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
            if ((length(method)!=5) || !is.numeric(method[3]) || !is.numeric(method[4])){
                stop("A subsetting condition for numeric-type data is not formatted properly.")
            }

            if (method[2]=="above") {
                left <- target_data > method[3]
            } else {
                left <- target_data >= method[3]
            }

            if (method[4]=="below") {
                right <- target_data < method[5]
            } else {
                right <- target_data <= method[5]
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
    # Then we can pass everything along to the normal function

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

##########

scdata <- readRDS(input_path("scdata"))
plotting_options <- input_json("plotting_options")
plot_setup <- input_json("plot_setup")

viz_fxn <- plot_setup$plot_type
plot_setup$plot_type <- NULL

# Parse reordering (Before input re-naming because of color.by -> var or color.var!)
if ("color_order" %in% names(plot_setup)) {

    if (length(plot_setup$order)>1) {
        # color_by is required to be discrete data to get here, which means it is a metadata
        # The method to set color orders for ggplot is via factor levels order
        scdata[[plot_setup$color_by]] <- factor(
            meta(plot_setup$color_by, scdata),
            levels = plot_setup$color_order,
        )
        plot_setup$order <- "increasing"
    } else {
        plot_setup$order <- plot_setup$color_order
    }

    plot_setup$color_order <- NULL
}

## Correct input names
names(plot_setup) <- gsub("_", ".", names(plot_setup))
rename_if_there <- function(current, proper) {
    if (current %in% names(plot_setup)) {
        ind <- which(names(plot_setup)==current)
        names(plot_setup)[ind] <- proper
    }
    plot_setup
}
plot_setup <- rename_if_there("x.by", "x.var")
plot_setup <- rename_if_there("y.by", "y.var")
if (viz_fxn=="dittoDimPlot") plot_setup <- rename_if_there("color.by", "var")
if (viz_fxn=="dittoScatterPlot") plot_setup <- rename_if_there("color.by", "color.var")
plot_setup <- rename_if_there("plot.title", "main")
plot_setup <- rename_if_there("plot.subtitle", "sub")
# BarPlot scaling (development choice diffs...)
plot_setup <- rename_if_there("scale.by", "scale")
if (!is.null(plot_setup$scale)) {
    plot_setup$scale <- switch(
        plot_setup$scale,
        NULL,
        'fraction' = 'percent',
        'counts' = 'count')
}

# # Convert from array/vector formats to singular value needed
# use_last <- function(elements) {
#     for (i in elements) {
#         if (i %in% names(plot_setup)) {
#             vals <- plot_setup[[i]]
#             plot_setup[[i]] <- vals[length(vals)]
#         }
#     }
#     plot_setup
# }
# plot_setup <- use_last(c("var", "y.var", "x.var", "color.by", "group.by"))

# Parse reduction_setup
if ("reduction.setup" %in% names(plot_setup)) {
    setup <- plot_setup$reduction.setup
    plot_setup$reduction.use <- setup[1]
    plot_setup$dim.1 <- as.numeric(setup[2])
    plot_setup$dim.2 <- as.numeric(setup[3])
    plot_setup$reduction.setup <- NULL

    # Replace recommendation
    if (plot_setup$reduction.use == "_Recommended_") {
        plot_setup$reduction.use <- plotting_options$Recommended_Reduction
    }
}

# Parse subsetting
if ( !is.null(plot_setup$cells.use) && length(plot_setup$cells.use)>0 ) {
    plot_setup$cells.use <- subset__scObj(scdata, plot_setup$cells.use)
} else {
    plot_setup$cells.use <- NULL
}

# Remove anything left as 'make' (so get filled with dittoSeq defaults!)
for (var in names(plot_setup)) {
    if (identical(plot_setup[[var]],"make")) {
        plot_setup[[var]] <- NULL
    }
}

# Add the dataset
plot_setup$object <- scdata

### Output the plot
# If do.hover=TRUE, turn that off for making a thumbnail
# Otherwise:
#   can make the thumbnail from the plot directly
#   split off legend so can ensure it won't override the main plot
if (!is.null(plot_setup$do.hover) && plot_setup$do.hover) {
    # Main plot with plotly
    # ToDo: Leave this up to the user and expose hover.data input!
    if (viz_fxn!="dittoBarPlot") {
        plot_setup$hover.data <- "var"
    }
    fig <- do.call(viz_fxn, plot_setup)
    fig_json <- plotly::plotly_json(fig, jsonedit = FALSE, pretty = FALSE)
    output_var(fig_json, 'plot.out')
    # Thumbnail
    plot_setup$do.hover <- FALSE
    ggsave(
        filename = output_path('plot.png'),
        plot = do.call(viz_fxn, plot_setup) + 
            theme_void() + theme(legend.position = "none"),
        units = "px",
        width = 300,
        height = 200
    )
} else {
    fig <- do.call(viz_fxn, plot_setup)
    ggsave(
        filename = output_path('plot.png'),
        plot = fig + 
            theme_void() + theme(legend.position = "none"),
        units = "px",
        width = 300,
        height = 200
    )
    ggsave(
        filename = output_path('plot.out'),
        # plot = dittoSeq:::.remove_legend(fig),
        plot = fig,
        device = "png",
        units = "px",
        width = 1200,
        height = 1000
    )
    # ggsave(
    #     filename = output_path('legend.png'),
    #     plot = dittoSeq:::.grab_legend(fig),
    #     device = "png",
    #     units = "in",
    #     width = 3,
    #     height = 6
    # )
}
