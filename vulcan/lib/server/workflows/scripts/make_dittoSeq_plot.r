library(Seurat)
library(dataflow)
library(dittoSeq)
library(plotly)

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
if (viz_fxn=="dittoDimPlot") {
    plot_setup <- rename_if_there("color.by", "var")
}
plot_setup <- rename_if_there("color.by", "color.var")
plot_setup <- rename_if_there("plot.title", "main")

# Convert from array/vector formats to singular value needed
use_last <- function(elements) {
    for (i in elements) {
        if (i %in% names(plot_setup)) {
            vals <- plot_setup[[i]]
            plot_setup[[i]] <- vals[length(vals)]
        }
    }
    plot_setup
}
plot_setup <- use_last(c("var", "y.var", "x.var", "color.by", "group.by"))

# Parse reduction_setup
if ("reduction.setup" %in% names(plot_setup)) {
    setup <- plot_setup$reduction.setup
    plot_setup$reduction.use <- setup[1]
    plot_setup$dim.1 <- as.numeric(setup[2])
    plot_setup$dim.2 <- as.numeric(setup[3])
    plot_setup$reduction.setup <- NULL
}

# Replace recommendations
if (plot_setup$reduction.use == "_Recommended_") {
    plot_setup$reduction.use <- plotting_options$Recommended_Reduction
}

# # Parse subsetting
# if (plot_setup$cells_use != {}) {
#     plot_setup$cells_use = subsetDF_index_targets(scdata, plot_setup$cells_use)
# } else {
    plot_setup$cells_use <- NULL
# }

# Add the dataset
plot_setup$object <- scdata

### Output png thumbnail
ggsave(
    filename = output_path('plot.png'),
    plot = do.call(viz_fxn, plot_setup) + 
        theme_void() + theme(legend.position = "none"),
    units = "px",
    width = 300,
    height = 200
)

# Turn on hover
# ToDo: Leave this up to the user and expose hover.data input!
plot_setup$do.hover <- TRUE
if (viz_fxn!="dittoBarPlot") {
    plot_setup$hover.data <- "var"
}

### Output plotly image
# Make & output plot
fig <- do.call(viz_fxn, plot_setup)
fig_json <- plotly::plotly_json(fig, jsonedit = FALSE, pretty = TRUE)
output_json(fig_json, 'plot.json')