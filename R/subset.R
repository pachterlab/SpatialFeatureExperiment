#' Subsetting SpatialFeatureExperiment objects
#'
#' The method for SFE reconstructs the spatial graphs when the SFE object is
#' subsetted as the \code{listw} objects encodes the nodes with indices which
#' are no longer valid after subsetting as some nodes are no longer present.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param i Row indices for subsetting.
#' @param j column indices for subsetting.
#' @param drop Logical. If \code{FALSE}, then a warning will be issued that the
#'   node indices in the graphs are no longer valid so the row and col graphs
#'   affected by subsetting are dropped. At present, this only works with the
#'   wrapper functions in this package that take in SFE objects and records the
#'   info required to reconstruct the graphs. While this argument is ignored for
#'   \code{SummarizedExperiment}
#' @param ... Passed to the \code{SingleCellExperiment} method of \code{[}.
#' @importFrom methods callNextMethod
#' @importFrom utils getFromNamespace
#' @return A subsetted \code{SpatialFeatureExperiment} object.
#' @name SpatialFeatureExperiment-subset
#' @aliases [,SpatialFeatureExperiment,ANY,ANY,ANY-method
#' @concept Non-spatial operations
#' @export
#' @examples
#' # Just like subsetting matrices and SingleCellExperiment
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe_subset <- sfe[seq_len(10), seq_len(10), drop = TRUE]
#' # Gives warning as graph reconstruction fails
#' \donttest{
#' sfe_subset <- sfe[seq_len(10), seq_len(10)]
#' }
setMethod(
    "[", c("SpatialFeatureExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        # Because the extra graphs and sample_ids result into invalid object
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }
        cn <- colnames(x)
        sample_ids0 <- sampleIDs(x)
        old_sample_colnames <- lapply(sample_ids0, function(s) {
            cn[colData(x)$sample_id %in% s]
        })
        # Subset the matrix and col and rowData
        # Suppress because `drop` will be used
        suppressWarnings(x <- callNextMethod())
        sample_ids <- sampleIDs(x)
        # Subset annotGeometries based on sample_id
        if (!is.null(annotGeometries(x))) {
            ag_sub <- annotGeometries(x)
            for (g in seq_along(ag_sub)) {
                ag_ind <- ag_sub[[g]]
                ag_ind <- ag_ind[ag_ind$sample_id %in% sample_ids, ]
                if (nrow(ag_ind) == 0L) ag_ind <- NULL
                ag_sub[[g]] <- ag_ind
            }
            annotGeometries(x) <- ag_sub
        }
        # Crop images with new bbox
        x <- .crop_imgs(x, bbox(x, sample_id = "all"))
        # Subset *Graphs based on sample_id and reconstruct row and colGraphs
        if (!is.null(spatialGraphs(x))) {
            graphs_sub <- int_metadata(x)$spatialGraphs
            graphs_sub <- graphs_sub[, names(graphs_sub) %in% sampleIDs(x),
                drop = FALSE
            ]
            if (!drop) {
                # Check which graphs need to be reconstructed
                # Wouldn't need reconstruction if the barcodes within one sample
                # are still in the same order
                cn2 <- colnames(x)
                new_sample_colnames <- lapply(sample_ids, function(s) {
                    cn2[colData(x)$sample_id %in% s]
                })
                old_sample_compare <- old_sample_colnames[sample_ids0 %in% sample_ids]
                samples_reconstruct <- mapply(
                    function(old, new) !isTRUE(all.equal(old, new)),
                    old = old_sample_compare,
                    new = new_sample_colnames,
                    SIMPLIFY = TRUE
                )
                for (s in which(samples_reconstruct)) {
                    for (m in seq_len(2)) { # Not reconstructing annotGraphs
                        # Not sure what to do differently with rowGraphs yet
                        for (g in seq_along(graphs_sub[[s]][[m]])) {
                            method_info <- attr(graphs_sub[[s]][[m]][[g]], "method")
                            if (is.null(method_info)) {
                                warning(
                                    "Graph reconstruction info is missing for sample ",
                                    names(graphs_sub)[s], " ", .margin_name(m), "Graph ",
                                    names(graphs_sub[[s]][[m]])[g], ". ",
                                    "Dropping graph.\n"
                                )
                                graphs_sub[[s]][[m]][[g]] <- NULL
                            } else {
                                if (requireNamespace(method_info$package[[1]], quietly = TRUE)) {
                                    fun <- getFromNamespace(method_info$FUN, method_info$package[[1]])
                                    if ("row.names" %in% names(method_info$args)) {
                                        method_info$args[["row.names"]] <-
                                            method_info$args[["row.names"]][j]
                                    }
                                    tryCatch(graphs_sub[[s]][[m]][[g]] <-
                                        do.call(fun, c(list(x = x), method_info$args)),
                                    error = function(e) {
                                        warning(
                                            "Graph reconstruction failed for sample ",
                                            names(graphs_sub)[s], " ",
                                            .margin_name(m), "Graph ",
                                            names(graphs_sub[[s]][[m]])[g],
                                            ": ", e, "Dropping graph.\n"
                                        )
                                        graphs_sub[[s]][[m]][[g]] <- NULL
                                    }
                                    )
                                } else {
                                    warning(
                                        "Package ", method_info$package[[1]],
                                        " used to construct graph for sample ",
                                        names(graphs_sub)[s], " ", .margin_name(m),
                                        "Graph ", names(graphs_sub[[s]][[m]])[g],
                                        " is not installed. ", "Dropping graph.\n"
                                    )
                                    graphs_sub[[s]][[m]][[g]] <- NULL
                                }
                            }
                        }
                    }
                }
                spatialGraphs(x) <- graphs_sub
            } else {
                message(
                    "Node indices in the graphs are no longer valid after subsetting. ",
                    "Dropping all row and col graphs."
                )
                spatialGraphs(x) <- graphs_sub
                spatialGraphs(x, MARGIN = 1) <- NULL
                spatialGraphs(x, MARGIN = 2) <- NULL
            }
        }
        validObject(x)
        return(x)
    }
)
