for (toolset in list.dirs("archimedes-r/libraries", full.names = TRUE, recursive = FALSE)) {
    install.packages(toolset, repos = NULL, type = "source")
}

