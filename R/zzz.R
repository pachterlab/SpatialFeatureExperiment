.onLoad <- function(libname, pkgname) {
    op <- options()
    if (!"SFE_graph_subset" %in% names(op)) {
        options(SFE_graph_subset = TRUE)
    }
    if (!"SFE_subset_crop" %in% names(op)) {
        options(SFE_subset_crop = TRUE)
    }
    if (!"SFE_subset_crop_max" %in% names(op)) {
        options(SFE_subset_crop_max = "100 MB")
    }
    invisible()
}
