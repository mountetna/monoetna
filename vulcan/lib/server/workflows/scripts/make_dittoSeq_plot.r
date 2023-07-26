dataflow::load_packages("dataflow", "dataSupport", "Seurat", "dittoSeq")

### This script...

scdata <- readRDS(input_path("scdata"))
plotting_options <- input_json("plotting_options")
plot_setup <- input_json("plot_setup")

viz_fxn <- plot_setup$plot_type
plot_setup$plot_type <- NULL

### Remove anything left as 'make' (so get filled with dittoSeq defaults!)
for (var in names(plot_setup)) {
    if (identical(plot_setup[[var]],"make") || (var=="vars_use" && identical(plot_setup[[var]],list())) ) {
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

## Inputs that don't map exactly
# 'split_adjust_free_y' actually translates to giving 'split.adjust = list(scale = "free_y")'
if ("split.adjust.free.y" %in% names(plot_setup)) {
    if (plot_setup$split.adjust.free.y) {
        plot_setup$split.adjust <- list(scale = "free_y")
    }
    plot_setup$split.adjust.free.y <- NULL
}

## Additional inputs always used
if (viz_fxn=="dittoBarPlot") {
    plot_setup$retain.factor.levels <- TRUE
}

### Add assay inputs for gene data
getAssay <- function(targ, object) {
    # Return: String, the assay name containing this gene
    for (assay_check in Seurat::Assays(object)) {
        if (isGene(targ, object, assay = assay_check)) {
            return(assay_check)
        }
    }
    stop(paste0("No assay found for ", targ))
}
addAssayInputIfGene <- function(gene_input_name, assay_input_name, setup = plot_setup, object = scdata) {
    # Return plot_setup with needed assay added to plot_setup[[assay_input_name]], when required
    if (gene_input_name %in% names(setup)) {
        gene_targ <- setup[[gene_input_name]]
        if (is(object, "Seurat") && !isMeta(gene_targ, object)) {
            setup[[assay_input_name]] <- getAssay(gene_targ, object)
        }
    }
    return(setup)
}
plot_setup <- addAssayInputIfGene("x.var", "assay.x")
plot_setup <- addAssayInputIfGene("y.var", "assay.y")
plot_setup <- addAssayInputIfGene("color.var", "assay.color")
plot_setup <- addAssayInputIfGene("color.var", "assay.hover")
plot_setup <- addAssayInputIfGene("var", "assay")
if (viz_fxn=="dittoDimPlot") plot_setup <- addAssayInputIfGene("var", "hover.assay")

### Parse reduction_setup
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

### Parse subsetting
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

### Turn off do_hover if too many points, using 100k as limit
if (viz_fxn %in% c('dittoDimPlot', 'dittoScatterPlot', 'dittoPlot') && ncol(scdata) > 25000)  {
    plot_setup$do.hover <- FALSE
}

### Add the dataset
plot_setup$object <- scdata

### Output the plot
# If do.hover=TRUE, turn that off for making a thumbnail
# Otherwise:
#   can make the thumbnail from the plot directly
#   ToDo: maybe split off legend so can ensure it won't override the main plot. Requires dev on the vulcan ui side to show these.
if (identical(plot_setup$do.hover, TRUE)) {
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
        dpi = 150,
        width = 1500,
        height = 1200
    )
    # Remove the Seurat 'object' from the plot's internals to keep it smaller
    fig$plot_env$object <- NULL
}

# Output the plot as file for editing
saveRDS(fig, file = output_path('plot.Rds'))
