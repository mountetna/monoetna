# For avoiding the linter's issue with `suppressPackageSatartupMessages()`
# For dev purposes, also adding a toggle to not suppress them
load_packages <- function(..., verbose = FALSE) {
    for (i in seq_len(...length())) {
        lib <- ...elt(i)
        if (verbose) {
            library(lib, character.only = TRUE)
        } else {
            suppressPackageStartupMessages({library(lib, character.only = TRUE)})
        }
    }

}