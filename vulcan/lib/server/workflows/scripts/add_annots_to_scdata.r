dataflow::load_packages("dataflow", "dataSupport", "Seurat", "dittoSeq")

scdata <- readRDS(input_path("scdata"))

### Pull in annotations by mapping to cells of their target clusters
annots <- read.csv(input_path("annots.csv"), header = TRUE)
clust_meta <- names(annots)[1]
inds_map <- match(scdata[[clust_meta, drop = TRUE]], annots[,1])
for (col in grep("^annot", names(annots)[-1], value = TRUE)) {
    scdata[[col]] <- annots[inds_map, col]
}

saveRDS(scdata, output_path("scdata"), compress = TRUE)