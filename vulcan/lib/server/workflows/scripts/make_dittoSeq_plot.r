dataflow::load_packages("dataflow", "dataSupport", "Seurat", "dittoSeq")

### This script...

scdata <- readRDS(input_path("scdata"))
plotting_options <- input_json("plotting_options")
plot_setup <- input_json("plot_setup")

viz_fxn <- plot_setup$plot_type
plot_setup$plot_type <- NULL

### Remove anything left as 'make' (so get filled with dittoSeq defaults!)
for (var in names(plot_setup)) {
    if (identical(plot_setup[[var]],"make")) {
        plot_setup[[var]] <- NULL
    }
}

### Replace "_Recommended_" metadata
## Alternatively, could have the UI serve these replacements, but I think that becomes visually finicky
replace_meta_rec <- function(value, meta_recommendations = plotting_options$Recommended_Metadata) {
    if (value %in% names(meta_recommendations)) {
        value <- meta_recommendations[[value]]
    }
    value
}
special_treatment <- c("cells_use")
for (i in seq_along(plot_setup)) {
    if (names(plot_setup)[i] %in% special_treatment) next

    if (length(plot_setup[[i]])>1) {
        for (j in seq_along(plot_setup[[i]])) {
            plot_setup[[i]][j] <- replace_meta_rec(plot_setup[[i]][j])
        }
    } else {
    plot_setup[i] <- replace_meta_rec(plot_setup[[i]])
    }
}

### Parse reorderings (Before input re-naming because of color.by -> either var or color.var!)
if ("color_order" %in% names(plot_setup)) {

    if (length(plot_setup$color_order)>1 || !plot_setup$color_order %in% c("increasing", "decreasing", "unordered", "randomize")) {
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
# (var_order & group_order are extra to dittoSeq features so require modifying the sc object itself!)
if ("var_order" %in% names(plot_setup)) {
    scdata@meta.data[,plot_setup$var] <- factor(
        scdata@meta.data[,plot_setup$var],
        levels = plot_setup$var_order
    )
    plot_setup$var_order <- NULL
}
if ("group_order" %in% names(plot_setup)) {
    scdata@meta.data[,plot_setup$group_by] <- factor(
        scdata@meta.data[,plot_setup$group_by],
        levels = plot_setup$group_order
    )
    plot_setup$group_order <- NULL
}

### Correct input names
# from '.'s to '_'s
names(plot_setup) <- gsub("_", ".", names(plot_setup))
rename_if_there <- function(current, proper) {
    if (current %in% names(plot_setup)) {
        ind <- which(names(plot_setup)==current)
        names(plot_setup)[ind] <- proper
    }
    plot_setup
}
# Other specifically different variable names
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

## Additional inputs always used
if (viz_fxn=="dittoBarPlot") {
    plot_setup$retain.factor.levels <- TRUE
}

# # Todo: ## Determine assays for gene data
# get_assay <- function(scdata, input, input_name, plot_setup = plot_setup) {
#     if (input_name %in% names(plot_setup) && !isMeta(plot_setup[[input_name]])) {
#
#     }
# }

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
    # Special-case recommendation replacement
    for (i in seq_along(plot_setup$cells.use$methods)) {
        method <- plot_setup$cells.use$methods[[i]]
        method[1] <- replace_meta_rec(method[1])
        plot_setup$cells.use$methods[[i]] <- method
    }
    # Parse
    plot_setup$cells.use <- subset__scObj(scdata, plot_setup$cells.use)
} else {
    plot_setup$cells.use <- NULL
}

# Add the dataset
plot_setup$object <- scdata

### Output the plot
# If do.hover=TRUE, turn that off for making a thumbnail
# Otherwise:
#   can make the thumbnail from the plot directly
#   ToDo: maybe split off legend so can ensure it won't override the main plot. Requires dev on the vulcan ui side to show these.
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
            theme_void() + theme(legend.position = "none") + ggtitle(NULL, NULL),
        units = "px",
        dpi = 50,
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
