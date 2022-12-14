.get_centroids <- function(x, type, MARGIN, sample_id) {
    if (type == "spatialCoords") {
        coords <- spatialCoords(x)[colData(x)$sample_id %in% sample_id, ]
        colnames(coords) <- c("x", "y")
        coords <- st_geometry(df2sf(coords))
    } else {
        # What to do with empty geometries?
        # Throw error for empty dimGeometries, since each item must have a geometry
        # Need a more helpful error message
        # Remove empty geometries for annotGeometries
        if (MARGIN < 3) {
            g <- dimGeometry(x, type, MARGIN, sample_id)
        } else {
            g <- annotGeometry(x, type, sample_id)
        }
        g <- .rm_empty_geometries(g, MARGIN)
        if (st_geometry_type(g, FALSE) == "POINT") {
            coords <- st_geometry(g)
        } else {
            coords <- st_centroid(st_geometry(g))
        }
        return(coords)
    }
}
.nb2listwdist2 <- function(nb, ...) UseMethod(".nb2listwdist2")

.nb2listwdist2.nbknn <- function(nb, type = "idw", style = "W", alpha = 1,
                                 dmax = NULL) {
    # Code adapted from spdep: https://github.com/r-spatial/spdep/blob/49c202d561da9565b0b70cf7462b7147feff59c2/R/nb2listwdist.R#L1
    distance <- attr(nb, "distance")
    if (type == "idw") {
        dmat <- distance^-alpha
        ind_finite <- is.finite(dmat)
        if (all(!ind_finite)) stop("All edge weights are infinite")
        if (any(!ind_finite)) {
            dmat[!ind_finite] <- max(dmat[ind_finite])
        }
    } else if (type == "exp") {
        dmat <- exp(-alpha * distance)
    } else if (type == "dpd") {
        dmat <- (1 - (distance/dmax)^alpha)^alpha
        dmat[dmat < 0] <- 0
    }

    if (!is.null(dmax) && dmax > 0) {
        dmat[distance > dmax] <- 0
    }

    # Row normalize the weights
    n <- nrow(dmat)
    if (style == "W") {
        dmat <- sweep(dmat, 1, rowSums(dmat), FUN = "/")
    } else if (style == "B") {
        dmat <- dmat > 0
    } else if (style %in% c("C", "U")) {
        D <- sum(dmat)
        if (style == "C") dmat <- n/D * dmat
        else dmat <- dmat/D
    } else if (style == "S") {
        q <- sqrt(rowSums(dmat^2))
        dmat <- sweep(dmat, 1, STATS = q, FUN = "/")
        Q <- sum(dmat)
        if (is.na(Q) || !(Q > 0))
            stop(paste("Failure in sum of intermediate weights:", Q))
        dmat <- dmat * n/Q
    } else if (style == "minmax") {
        mm <- min(max(rowSums(dmat)), max(colSums(dmat)))
        dmat <- dmat / mm
    }

    if (anyNA(dmat)) {
        stop("NAs in coding scheme weights list")
    }

    # Construct the listw, should decouple finding edges and weights in future version
    glist <- asplit(dmat, 1)
    glist <- lapply(seq_along(glist), function(i) glist[[i]][ord[[i]]])

    listw <- list(style = style,
                  neighbours = nb,
                  weights = glist)
    class(listw) <- "listw"
    attr(listw, "region.id") <- attr(nb, "region.id")
    listw
}

#' @importFrom spdep card
.nb2listwdist2.nbdnn <- function(nb, type = "idw", style = "W", alpha = 1,
                                 dmax = NULL) {
    # Code adapted from spdep: https://github.com/r-spatial/spdep/blob/49c202d561da9565b0b70cf7462b7147feff59c2/R/nb2listwdist.R#L1
    cardnb <- card(nb)
    vlist <- vector("list", length = length(nn$index))
    glist <- attr(nb, "distance")
    if (type == "idw") {
        for (i in 1:n) {
            if (cardnb[i] > 0) {
                vlist[[i]] <- glist[[i]]^-alpha
                if(!is.null(dmax) && dmax > 0)
                    vlist[[i]][which(glist[[i]] > dmax)] <- 0
            }
        }
        uvlist <- unlist(vlist)
        fins <- is.finite(uvlist)
        if (all(!fins)) stop("no finite general weights")
        if (any(!fins)) {
            max_finite <- max(uvlist[fins])
            for(i in 1:n) {
                vlist[[i]][is.infinite(vlist[[i]])] <- max_finite
            }
        }
    }

    if (type == "exp") {
        for (i in 1:n) {
            if (cardnb[i] > 0) {
                vlist[[i]] <- exp(glist[[i]] * ((-1) * alpha))
                if(!is.null(dmax) && dmax > 0)
                    vlist[[i]][glist[[i]] > dmax] <- 0
            }
        }
    }

    if (type == "dpd") {
        if (is.null(dmax)) stop("DPD weights require a maximum distance threshold")
        if (dmax <= 0) stop("DPD weights require a positive maximum distance threshold")
        for (i in 1:n) {
            if (cardnb[i] > 0) {
                vlist[[i]] <- (1 - (glist[[i]] / dmax)^alpha)^alpha
                vlist[[i]][vlist[[i]] < 0] <- 0
            }
        }
    }

    if(style != "raw") {

        if (zero.policy) {
            eff.n <- n - sum(cardnb == 0)
            if (eff.n < 1) stop("No valid observations")
        } else eff.n <- n

        if (style == "W") {
            d <- unlist(lapply(glist, sum))
            for (i in 1:n) {
                if (cardnb[i] > 0) {
                    if (d[i] > 0) vlist[[i]] <- glist[[i]] / d[i]
                    else vlist[[i]] <- 0
                }
            }
            attr(vlist, "comp") <- list(d=d)
        }

        if (style == "B") {
            for (i in 1:n) {
                if (cardnb[i] > 0) vlist[[i]] <- as.numeric(glist[[i]] > 0)
            }
        }

        if (style == "C" || style == "U") {
            D <- sum(unlist(glist))
            if (is.na(D) || !(D > 0))
                stop(paste("Failure in sum of weights:", D))
            for (i in 1:n) {
                if (cardnb[i] > 0) {
                    if (style == "C")
                        vlist[[i]] <- (eff.n/D) * glist[[i]]
                    else
                        vlist[[i]] <- (1/D) * glist[[i]]
                }
            }
        }

        if (style == "S") {
            glist2 <- lapply(glist, function(x) x^2)
            q <- sqrt(unlist(lapply(glist2, sum)))
            for (i in 1:n) {
                if (cardnb[i] > 0) {
                    if (q[i] > 0) glist[[i]] <- (1/q[i]) * glist[[i]]
                    else glist[[i]] <- 0
                }
            }
            Q <- sum(unlist(glist))
            if (is.na(Q) || !(Q > 0))
                stop(paste("Failure in sum of intermediate weights:", Q))
            for (i in 1:n) {
                if (cardnb[i] > 0)
                    vlist[[i]] <- (eff.n/Q) * glist[[i]]
            }
            attr(vlist, "comp") <- list(q=q, Q=Q, eff.n=eff.n)
        }
    }

    if (!zero.policy)
        if (any(is.na(unlist(vlist))))
            stop ("NAs in coding scheme weights list")

    if (style == "minmax") {
        res <- list(style=style, neighbours=neighbours, weights=vlist)
        class(res) <- c("listw", "nb")
        mm <- minmax.listw(res)
        vlist <- lapply(vlist, function(x) (1/c(mm)) * x)
    }

    listw <- list(style = style,
                  neighbours = nb,
                  weights = vlist)
    class(listw) <- "listw"
    attr(listw, "region.id") <- attr(nb, "region.id")
    listw
}

#' @importFrom BiocNeighbors AnnoyParam KmknnParam
.knn_bioc <- function(coords, k = 1, BNPARAM = AnnoyParam(),
                      BPPARAM = SerialParam(), row.names = NULL) {
    nn <- findKNN(coords, k = k, BNPARAM = BNPARAM, BPPARAM = BPPARAM)
    # Split by row
    nb <- asplit(nn$index, 1)

    # Sort based on index
    ord <- lapply(nb, order)
    nb <- lapply(seq_along(nb), function(i) nb[[i]][ord[[i]]])
    attr(nb, "distance") <- nn$distance
    attr(nb, "region.id") <- row.names
    class(nb) <- c("nb", "nbknn")
    nb
}

.dnn_bioc <- function(coords, d2, BNPARAM = KmknnParam(),
                      BPPARAM = SerialParam(), row.names = NULL,
                      zero.policy = TRUE) {
    nn <- findNeighbors(coords, threshold = d2, BNPARAM = BNPARAM, BPPARAM = BPPARAM)

    # I might pad with 0 and convert the whole thing into a matrix
    # and use the same matrix based code in .knn_bioc
    # If that way is faster. But I'm not optimizing yet.
    cardnb <- lengths(nn$index)
    if (!zero.policy)
        if (any(cardnb == 1)) stop("Empty neighbour sets found")

    # Remove self from index and reorder
    nb <- glist <- vector("list", length = length(nn$index))
    n <- length(nb)
    for (i in seq_along(nb)) {
        index <- nn$index[[i]]
        ind_use <- index != i
        if (any(ind_use)) {
            v <- index[ind_use]
            ord <- order(v)
            nb[[i]] <- v[ord]
            glist[[i]] <- nn$distance[[i]][ind_use][ord]
        } else {
            nb[[i]] <- 0L
            glist[[i]] <- 0L
        }
    }
    attr(nb, "distance") <- glist
    attr(nb, "region.id") <- row.names
    class(nb) <- c("nb", "nbdnn")
    nb
}

.knn_sfe <- function(coords, k = 1, row.names = NULL, nn_method = "bioc",
                     use_kd_tree = TRUE, BNPARAM = KmknnParam(),
                     BPPARAM = SerialParam()) {
    if (nn_method == "spdep")
        knn2nb(knearneigh(coords, k = k, use_kd_tree = use_kd_tree),
               row.names = row.names)
    else
        .knn_bioc(coords, k = k, BPPARAM = BPPARAM, BNPARAM = BNPARAM)
}

.dnn_sfe <- function(coords, d1, d2, row.names = NULL, nn_method = "bioc",
                     use_kd_tree = TRUE, BNPARAM = KmknnParam(),
                     BPPARAM = SerialParam(), zero.policy = TRUE) {
    if (nn_method == "spdep")
        dnearneigh(coords, d1, d2, use_kd_tree = use_kd_tree, row.names = row.names)
    else
        .dnn_bioc(coords, d2 = d2, BNPARAM = BNPARAM, BPPARAM = BPPARAM,
                  row.names = row.names, zero.policy = zero.policy)
}
.g2nb_sfe <- function(coords, fun, nnmult = 3, sym = FALSE, row.names = NULL) {
    # Either gabrielneigh or relativeneigh
    g <- fun(coords, nnmult)
    graph2nb(g, sym = sym, row.names = row.names)
}
.gabriel_sfe <- function(coords, nnmult = 3, sym = FALSE, row.names = NULL) {
    .g2nb_sfe(coords, gabrielneigh, nnmult, sym, row.names)
}
.relative_sfe <- function(coords, nnmult = 3, sym = FALSE, row.names = NULL) {
    .g2nb_sfe(coords, relativeneigh, nnmult, sym, row.names)
}
.soi_sfe <- function(coords, quadsegs = 10, sym = FALSE, row.names = NULL) {
    g <- soi.graph(tri2nb(coords), coords, quadsegs)
    graph2nb(g, sym = sym, row.names = row.names)
}

.comp_graph_sample <- function(x, sample_id, type, MARGIN, method, dist_type,
                               args, extra_args_use, glist, style, zero.policy,
                               alpha, dmax, fun_use) {
    if (!"row.names" %in% names(args) &&
        "row.names" %in% extra_args_use && MARGIN < 3) {
        args$row.names <- colnames(x)[colData(x)$sample_id == sample_id]
    }
    if (method != "poly2nb") {
        coords <- .get_centroids(x, type, MARGIN, sample_id)
        nb_out <- do.call(fun_use, c(list(coords = coords), args))
    } else {
        if (MARGIN < 3) {
            coords <- dimGeometry(x, type, MARGIN, sample_id)
        } else {
            coords <- annotGeometry(x, type, sample_id)
        }
        coords <- .rm_empty_geometries(coords, MARGIN)
        if (st_geometry_type(coords, FALSE) != "POLYGON") {
            stop("poly2nb can only be used on POLYGON geometries.")
        }
        nb_out <- do.call(fun_use, c(list(pl = coords), args))
    }
    if (dist_type == "none") {
        if (style == "raw") style <- "W"
        out <- nb2listw(nb_out, glist, style, zero.policy)
    } else {
        out <- nb2listwdist(nb_out, coords,
            type = dist_type, style = style,
            longlat = FALSE, zero.policy = zero.policy,
            alpha = alpha, dmax = dmax
        )
    }

    attr(out, "method") <- list(
        FUN = "findSpatialNeighbors",
        package = "SpatialFeatureExperiment",
        args = c(
            method = method, args,
            dist_type = dist_type,
            glist = glist,
            style = style,
            alpha = alpha,
            dmax = dmax,
            zero.policy = zero.policy,
            sample_id = sample_id,
            type = type,
            MARGIN = MARGIN
        )
    )
    out
}

#' Find spatial neighborhood graph
#'
#' This function wraps all spatial neighborhood graphs implemented in the
#' package \code{spdep} for the \code{SpatialFeatureExperiment} (SFE) class, to
#' find spatial neighborhood graphs for the entities represented by columns or
#' rows of the gene count matrix in the SFE object or spatial entities in the
#' \code{annotGeometries} field of the SFE object. Results are stored as
#' \code{listw} objects in the \code{spatialGraphs} field of the SFE object, as
#' \code{listw} is used in many methods that facilitate the spatial neighborhood
#' graph in the \code{spdep}, \code{spatialreg}, and \code{adespatial}. The edge
#' weights of the graph in the \code{listw} object are by default style W (see
#' \code{\link{nb2listw}}) and the unweighted neighbor list is in the
#' \code{neighbours} field of the \code{listw} object.
#'
#' @inheritParams spdep::nb2listw
#' @param x A \code{\link{SpatialFeatureExperiment}} object.
#' @param MARGIN Just like in \code{\link{apply}}, where 1 stands for row, 2
#'   stands for column. Here, in addition, 3 stands for annotation, to query the
#'   \code{\link{annotGeometries}}, such as nuclei segmentation in a Visium data
#' @param sample_id Which sample(s) in the SFE object to use for the graph. Can
#'   also be "all", which means this function will compute the graph for all
#'   samples independently.
#' @param type Name of the geometry associated with the MARGIN of interest for
#'   which to compute the graph.
#' @param method Name of function in the package \code{spdep} to use to find the
#'   spatial neighborhood graph.
#' @param dist_type Type of distance-based weight. "none" means not using
#'   distance-based weights; the edge weights of the spatial neighborhood graph
#'   will be entirely determined by the \code{style} argument. "idw" means
#'   inverse distance weighting. "exp" means exponential decay. "dpd" means
#'   double-power distance weights. See \code{\link[spdep]{nb2listwdist}} for
#'   details.
#' @param nn_method Method to find k nearest neighbors and distance based
#' neighbors. Can be either "bioc" or "spdep". For "bioc", methods from
#' \code{BiocNeighbors} are used. For "spdep", methods from the \code{spdep}
#' package are used. The "bioc" option is more scalable to larger datasets and
#' supports multithreading.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object for multithreading.
#' Only used for k nearest neighbor and distance based neighbor with
#' \code{nn_method = "bioc"}.
#' @param BNPARAM A \code{\link{BiocNeighborParam}} object specifying the
#' algorithm to find k nearest neighbors and distance based neighbors with
#' \code{nn_method = "bioc"}. For distance based neighbors, only
#' \code{\link{KmknnParam}} and \code{\link{VptreeParam}} are applicable.
#' @param alpha Only relevant when \code{dist_type = "dpd"}.
#' @param dmax Only relevant when \code{dist_type = "dpd"}.
#' @param ... Extra arguments passed to the \code{spdep} function stated in the
#'   \code{method} argument, such as \code{k}, \code{use_kd_tree}, \code{d1},
#'   \code{d2}, \code{nnmult}, \code{sym}, and \code{quadsegs}. Note that any
#'   arguments about using longitude and latitude, which are irrelevant, are
#'   ignored.
#' @return For one sample, then a \code{listw} object representing the graph,
#'   with an attribute "method" recording the function used to build the graph,
#'   its arguments, and information about the geometry for which the graph was
#'   built. The attribute is used to reconstruct the graphs when the SFE object
#'   is subsetted since some nodes in the graph will no longer be present. If
#'   sample_id = "all" or has length > 1, then a named list of \code{listw}
#'   objects, whose names are the sample_ids. To add the list for multiple
#'   samples to a SFE object, specify the \code{name} argument in the
#'   \code{\link{spatialGraphs}} replacement method, so graph of the same name
#'   will be added to the SFE object for each sample.
#' @importFrom spdep tri2nb knearneigh dnearneigh gabrielneigh relativeneigh
#'   soi.graph knn2nb graph2nb nb2listw poly2nb nb2listwdist
#' @aliases findSpatialNeighbors
#' @note \code{style = "raw"} is only applicable when \code{dist_type} is not
#'   "none". If \code{dist_type = "none"} and \code{style = "raw"}, then style
#'   will default to "W". Using distance based weights does not supplant finding
#'   a spatial neighborhood graph. The spatial neighborhood graph is first found
#'   and then its edges weighted based on distance in this function.
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' # sample_id is optional when only one sample is present
#' g <- findSpatialNeighbors(sfe, sample_id = "Vis5A")
#' attr(g, "method")
#' # Returns named list for multiple samples
#' sfe2 <- McKellarMuscleData(dataset = "small2")
#' sfe_combined <- cbind(sfe, sfe2)
#' gs <- findSpatialNeighbors(sfe, sample_id = "all")
setMethod(
    "findSpatialNeighbors", "SpatialFeatureExperiment",
    function(x, sample_id = NULL, type = "spatialCoords", MARGIN = 2,
             method = c(
                 "tri2nb", "knearneigh", "dnearneigh",
                 "gabrielneigh", "relativeneigh", "soi.graph",
                 "poly2nb"
             ),
             dist_type = c("none", "idw", "exp", "dpd"),
             glist = NULL, style = c(
                 "raw", "W", "B", "C", "U",
                 "minmax", "S"
             ), nn_method = c("bioc", "spdep"),
             alpha = 1, dmax = NULL, zero.policy = NULL, ...) {
        method <- match.arg(method)
        dist_type <- match.arg(dist_type)
        style <- match.arg(style)
        sample_id <- .check_sample_id(x, sample_id, one = FALSE)
        extra_args_use <- switch(method,
            tri2nb = "row.names",
            knearneigh = c("k", "use_kd_tree"),
            dnearneigh = c(
                "d1", "d2", "use_kd_tree",
                "row.names"
            ),
            gabrielneigh = c(
                "nnmult", "sym",
                "row.names"
            ),
            relativeneigh = c(
                "nnmult", "sym",
                "row.names"
            ),
            soi.graph = c(
                "quadsegs", "sym",
                "row.names"
            ),
            poly2nb = c(
                "row.names", "snap", "queen",
                "useC", "foundInBox"
            )
        )
        args <- list(...)
        args <- args[names(args) %in% extra_args_use]
        fun_use <- switch(method,
            tri2nb = tri2nb,
            knearneigh = .knn_sfe,
            dnearneigh = .dnn_sfe,
            gabrielneigh = .gabriel_sfe,
            relativeneigh = .relative_sfe,
            soi.graph = .soi_sfe,
            poly2nb = poly2nb
        )
        if (length(sample_id) == 1L) {
            out <- .comp_graph_sample(
                x, sample_id, type, MARGIN, method,
                dist_type, args, extra_args_use, glist,
                style, zero.policy, alpha, dmax, fun_use
            )
        } else {
            out <- lapply(sample_id, function(s) {
                .comp_graph_sample(
                    x, s, type, MARGIN, method, dist_type,
                    args, extra_args_use, glist, style,
                    zero.policy, alpha, dmax, fun_use
                )
            })
            names(out) <- sample_id
        }
        return(out)
    }
)

.comp_visium_graph <- function(x, sample_id, style, zero.policy) {
    bcs_use <- colnames(x)[colData(x)$sample_id == sample_id]
    bcs_use2 <- sub("[-\\d]+$", "", bcs_use, perl = TRUE)
    # visium_row_col <- SpatialFeatureExperiment::visium_row_col
    coords_use <- visium_row_col[
        match(bcs_use2, visium_row_col$barcode),
        c("col", "row")
    ]
    # So adjacent spots are equidistant
    coords_use$row <- coords_use$row * sqrt(3)
    g <- dnearneigh(as.matrix(coords_use),
        d1 = 1.9, d2 = 2.1,
        row.names = bcs_use
    )
    out <- nb2listw(g, style = style, zero.policy = zero.policy)
    attr(out, "method") <- list(
        FUN = "findVisiumGraph",
        package = "SpatialFeatureExperiment",
        args = list(
            style = style,
            zero.policy = zero.policy,
            sample_id = sample_id
        )
    )
    out
}

#' Find spatial neighborhood graphs for Visium spots
#'
#' Visium spots are arranged in a hexagonal grid. This function uses the known
#' locations of the Visium barcodes to construct a neighborhood graph, so
#' adjacent spots are connected by edges. Since the known rows and columns of
#' the spots are used, the unit the spot centroid coordinates are in does not
#' matter.
#'
#' @inheritParams spdep::nb2listw
#' @param x A \code{SpatialFeatureExperiment} object with Visium data. Column
#'   names of the gene count matrix must be Visium barcodes, which may have a
#'   numeric suffix to distinguish between samples (e.g. "AAACAACGAATAGTTC-1").
#' @param sample_id Which sample(s) in the SFE object to use for the graph. Can
#'   also be "all", which means this function will compute the graph for all
#'   samples independently.
#' @importFrom spdep dnearneigh nb2listw
#' @return For one sample, then a \code{listw} object representing the graph,
#'   with an attribute "method" recording the function used to build the graph,
#'   its arguments, and information about the geometry for which the graph was
#'   built. The attribute is used to reconstruct the graphs when the SFE object
#'   is subsetted since some nodes in the graph will no longer be present. If
#'   sample_id = "all" or has length > 1, then a named list of \code{listw}
#'   objects, whose names are the sample_ids. To add the list for multiple
#'   samples to a SFE object, specify the \code{name} argument in the
#'   \code{\link{spatialGraphs}} replacement method, so graph of the same name
#'   will be added to the SFE object for each sample.
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' g <- findVisiumGraph(sfe)
#' # For multiple samples, returns named list
#' sfe2 <- McKellarMuscleData(dataset = "small2")
#' sfe_combined <- cbind(sfe, sfe2)
#' gs <- findVisiumGraph(sfe, sample_id = "all")
findVisiumGraph <- function(x, sample_id = NULL, style = "W",
                            zero.policy = NULL) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (length(sample_id) == 1L) {
        out <- .comp_visium_graph(x, sample_id, style, zero.policy)
    } else {
        out <- lapply(
            sample_id,
            function(s) .comp_visium_graph(x, s, style, zero.policy)
        )
        names(out) <- sample_id
    }
    return(out)
}
