.onLoad <- function(libname, pkgname) {
    op <- options()
    if (!"SFE_graph_subset" %in% names(op)) {
        options(SFE_graph_subset = TRUE)
    }
    invisible()
}
