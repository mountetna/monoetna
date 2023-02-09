library(dataflow)
library(SeuratDisk)

tmp_dir <- tempdir()
h5ad_path <- file.path(tmp_dir, 'scdata.h5ad')
h5seurat <- file.path(tmp_dir, 'scdata.h5seurat')

file.copy(input_path("scdata.h5ad"), h5ad_path)

Convert(h5ad_path, dest = "h5seurat")

if (file.size(h5seurat)==0) {
    stop("scanpy to Seurat conversion error: SeuratDisk::Convert() yielded a 0 byte file.")
}

scdata <- LoadH5Seurat(h5seurat)

saveRDS(scdata, output_path("scdata.Rds"))