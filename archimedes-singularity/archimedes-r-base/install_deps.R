install.packages("remotes")

install.packages("BiocManager")
BiocManager::install(version = "3.19")
stopifnot(BiocManager::version()=='3.19')

cran_deps <- list(
    c("jsonlite", "1.8.8"),
    c("magmaR"),
    c("Seurat", "5.1.0"),
    c("plotly", "4.10.4")
)
for (set in cran_deps) {
    pkg <- set[1]
    if (length(set) > 1 && !is.null(set[2])) {
        ver <- set[2]
        remotes::install_version(pkg, version = ver, upgrade = FALSE)
    } else {
        install.packages(pkg)
    }
}

github_deps <- list(
    c("mojaveazure/seurat-disk", "877d4e18ab38c686f5db54f8cd290274ccdbe295"),
    c("dtm2451/dittoSeq"),
    c("dtm2451/dittoViz")
)
for (set in github_deps) {
    repo <- set[1]
    ref <- ifelse(length(set) > 1 && !is.null(set[2]), set[2], "HEAD")
    BiocManager::install(repo, ref = ref, update = FALSE)
}