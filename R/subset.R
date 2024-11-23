.is0 <- function(x) {
    if (length(x)) {
        if (is.logical(x)) !any(x) else FALSE
    } else TRUE
}

.subset_listw <- function(x, subset, zero.policy=attr(x, "zero.policy"), ...) {
    nb <- x$neighbours
    vlist <- x$weights
    style <- x$style
    m_wts <- as_dgRMatrix_listw(x)
    m_wts <- m_wts[subset, subset]
    mat2listw(m_wts, style = style, zero.policy = zero.policy)
}

#' Subsetting SpatialFeatureExperiment objects
#'
#' The SFE method has special treatment for the spatial graphs. In \code{listw},
#' the neighbors are indicated by indices, which will change after subsetting.
#' The \code{SFE_graph_subset} option determines whether the graphs are
#' subsetted or reconstructed. In the default (\code{options(SFE_graph_subset =
#' TRUE)}), the graphs are subsetted, in which case singletons may be produced.
#' For \code{options(SFE_graph_subset = FALSE)}, which is the behavior of
#' versions earlier than Bioc 3.20, the graphs are reconstructed with the
#' parameters recorded in an attribute of the graphs. This option can result
#' into different graphs. For example, suppose we start with a k nearest
#' neighbor graph. After subsetting, cells at the boundary of the region used to
#' subset the SFE object may lose some of their neighbors. In contrast, when the
#' graph is reconstructed, these same edge cells will gain other cells that
#' remain after subsetting as neighbors in the new KNN graph.
#'
#' The option \code{SFE_graph_subset} was introduced because subsetting is
#' usually faster than reconstructing and in some cases such as distance-based
#' neighbors and Visium spot adjacency give the same results. It was introduced
#' also because of the development of \code{alabster.sfe} for a
#' language-agnostic on-disk serialization of SFE objects and some parameters
#' used to construct graphs have special classes whose \code{alabaster} methods
#' have not been implemented, such as \code{BPPARAM} and \code{BNPARAM}, so when
#' reconstructing, the defaults for those arguments will be used.
#'
#' The edge weights will be recomputed from the binary neighborhood indicator
#' with the same normalization style as the original graph, such as "W" for row
#' normalization. When distance-based edge weights are used instead of the
#' binary indicator, the edge weights will be re-normalized, which is mostly
#' some rescaling. This should give the same results as recomputing the distance
#' based edge weights for styles "raw", "W", and "B" since the distances
#' themselves don't change, but the effects of other more complicated styles of
#' re-normalization on spatial statistics should be further investigated.
#'
#' @param x A \code{SpatialFeatureExperiment} object.
#' @param i Row indices for subsetting.
#' @param j column indices for subsetting.
#' @param drop Only used if graphs are reconstructed
#'   (\code{options(SFE_graph_subset = FALSE)}). If \code{TRUE} then
#'   \code{colGraphs} are dropped but \code{annotGraphs} are kept.
#' @param ... Passed to the \code{SingleCellExperiment} method of \code{[}.
#' @importFrom methods callNextMethod
#' @importFrom utils getFromNamespace
#' @importFrom spdep mat2listw
#' @importFrom spatialreg as_dgRMatrix_listw
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
        if (length(annotGeometries(x))) {
            ag_sub <- annotGeometries(x)
            for (g in seq_along(ag_sub)) {
                ag_ind <- ag_sub[[g]]
                ag_ind <- ag_ind[ag_ind$sample_id %in% sample_ids, ]
                ag_sub[[g]] <- ag_ind
            }
            annotGeometries(x) <- ag_sub
        }
        # Remove rowGeometries when an entire sample got removed
        if (!is.null(rowGeometryNames(x))) {
            samples_rm <- setdiff(sample_ids0, sample_ids)
            if (length(samples_rm)) {
                rowGeometries(x, sample_id = samples_rm) <- NULL
            }
        }
        # Crop images with new bbox
        if (!missing(j)) {
            x <- .crop_imgs(x, bbox(x, sample_id = "all"))
        }
        # Subset *Graphs based on sample_id and reconstruct row and colGraphs
        if (!is.null(spatialGraphs(x)) && (!missing(j) && !.is0(j))) {
            do_subset <- getOption("SFE_graph_subset")
            graphs_sub <- int_metadata(x)$spatialGraphs
            graphs_sub <- graphs_sub[, names(graphs_sub) %in% sampleIDs(x),
                drop = FALSE
            ]
            if (drop && !do_subset) {
                message(
                    "Node indices in the graphs are no longer valid after subsetting. ",
                    "Dropping all row and col graphs."
                )
                spatialGraphs(x) <- graphs_sub
                spatialGraphs(x, MARGIN = 1) <- NULL
                spatialGraphs(x, MARGIN = 2) <- NULL
                validObject(x)
                return(x)
            }
            # Check which graphs need to be subsetted
            # Wouldn't need reconstruction if the barcodes within one sample
            # are still in the same order
            cn2 <- colnames(x)
            new_sample_colnames <- lapply(sample_ids, function(s) {
                cn2[colData(x)$sample_id %in% s]
            })
            old_sample_compare <- old_sample_colnames[sample_ids0 %in% sample_ids]
            samples_subset <- mapply(
                function(old, new) !isTRUE(all.equal(old, new)),
                old = old_sample_compare,
                new = new_sample_colnames,
                SIMPLIFY = TRUE
            )
            for (s in which(samples_subset)) {
                j_sample <- old_sample_compare[[s]] %in% new_sample_colnames[[s]]
                if (sum(j_sample) == 1) {
                    message("Only 1 column left in sample ", sample_ids[[s]],
                            ", graphs are meaningless, dropping graphs.")
                    graphs_sub[[s]] <- list(row = NULL, col = NULL, annot = NULL)
                    next
                }
                for (m in seq_len(2)) { # Not subsetting annotGraphs
                    # Not sure what to do differently with rowGraphs yet
                    for (g in seq_along(graphs_sub[[s]][[m]])) {
                        method_info <- attr(graphs_sub[[s]][[m]][[g]], "method")
                        if (do_subset) {
                            gr <- .subset_listw(graphs_sub[[s]][[m]][[g]], j_sample)
                            attr(gr, "method") <- method_info
                            graphs_sub[[s]][[m]][[g]] <- gr
                        } else {
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
            }
            spatialGraphs(x) <- graphs_sub
        }
        if (!missing(j) && .is0(j)) spatialGraphs(x) <- NULL
        validObject(x)
        return(x)
    }
)
