#' Simple geometry predicates
#'
#' Unlike functions in \code{sf} like \code{st_intersects}, this function simply
#' returns a logical vector indicating whether each geometry in \code{x}
#' intersects (or returns \code{TRUE} from other predicates) anything in
#' \code{y}, preferably when \code{y} only contains a small number of geometries
#' or is one single MULTI geometry. This is useful when cropping or subsetting
#' an SFE object with a geometry, such as tissue boundary or histological region
#' polygons or a bounding box.
#'
#' @param x An object of class \code{sf}, \code{sfc}, or \code{sfg}.
#' @param y Another object of class \code{sf}, \code{sfc}, or \code{sfg}.
#' @param pred A geometric binary predicate function, such as
#'   \code{\link{st_intersects}}. It should return an object of class
#'   \code{sgbp}, for sparse predicates.
#' @return For \code{st_any_*}, a logical vector indicating whether each
#'   geometry in \code{x} intersects (or other predicates such as is covered by)
#'   anything in \code{y}. Simplified from the \code{sgbp} results which
#'   indicate which item in \code{y} each item in \code{x} intersects, which
#'   might not always be relevant. For \code{st_n_*}, an integer vector
#'   indicating the number of geometries in y returns TRUE for each geometry in
#'   x.
#' @concept Geometric operations
#' @export
#' @importFrom sf st_intersects st_agr<- st_drop_geometry st_as_sfc st_cast
#'   st_is_empty st_disjoint st_z_range st_zm
#' @importFrom stats aggregate
#' @examples
#' library(sf)
#' pts <- st_sfc(
#'     st_point(c(.5, .5)), st_point(c(1.5, 1.5)),
#'     st_point(c(2.5, 2.5))
#' )
#' pol <- st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
#' st_any_pred(pts, pol, pred = st_disjoint)
#' st_any_intersects(pts, pol)
#' st_n_pred(pts, pol, pred = st_disjoint)
#' st_n_intersects(pts, pol)
st_any_pred <- function(x, y, pred) lengths(pred(x, y)) > 0L
# TODO: put the item with more geometries in the x position. This way is much faster.
# If swapping positions, it shouldn't be hard to recover it if the pred is symmetric

#' @rdname st_any_pred
#' @export
st_any_intersects <- function(x, y) st_any_pred(x, y, st_intersects)

#' @rdname st_any_pred
#' @export
st_n_pred <- function(x, y, pred) lengths(pred(x, y))

#' @rdname st_any_pred
#' @export
st_n_intersects <- function(x, y) st_n_pred(x, y, st_intersects)

.crop_geometry <- function(g, y, samples_use, op, sample_col = NULL,
                           id_col = "id", remove_empty = FALSE) {
    if (!nrow(g)) return(g)
    # g may have column sample_id as in colGeometry and annotGeometry
    # g should also have row names if it is colGeometry
    if (!id_col %in% names(g)) {
        rm_id <- TRUE
        if (!is.null(rownames(g))) {
            g[[id_col]] <- rownames(g)
        } else {
            g[[id_col]] <- as.character(seq_len(nrow(g)))
        }
    } else {
        rm_id <- FALSE
    }
    if (!is.null(sample_col)) g$sample_id <- sample_col
    gs <- split(g, g$sample_id)
    gs_sub <- lapply(names(gs), function(s) {
        if (s %in% samples_use) {
            if ("sample_id" %in% names(y)) {
                y_use <- st_geometry(y[y$sample_id == s, ])
            } else {
                y_use <- st_geometry(y)
            }
            .g <- gs[[s]][, c("geometry", id_col)]
            st_agr(.g) <- "constant"
            if (!is.null(st_z_range(.g)))
                y_use <- st_zm(y_use, drop = FALSE, what = "Z")
            o <- op(.g, y_use)
            # Aggregate in case cropping broke some items into multiple pieces
            if (any(!rownames(o) %in% rownames(.g))) {
                o <- aggregate(o,
                               by = setNames(list(id = o[[id_col]]), id_col),
                               FUN = unique
                )
            }
            return(merge(o, st_drop_geometry(gs[[s]]), by = id_col, all = TRUE))
        } else {
            return(gs[[s]])
        }
    })
    gs_sub <- do.call(rbind, gs_sub)
    # Convert st_GEOMETRY to a more specific type
    if (st_geometry_type(gs_sub, by_geometry = FALSE) == "GEOMETRY") {
        gs_sub <- st_cast(gs_sub)
    }
    gs_sub <- gs_sub[, names(g)]
    rownames(gs_sub) <- rownames(g)[match(gs_sub[[id_col]], g[[id_col]])]
    if (rm_id) gs_sub[[id_col]] <- NULL
    if (remove_empty) gs_sub <- gs_sub[!st_is_empty(gs_sub), ]
    gs_sub
}

.annot_fun <- function(x, y, colGeometryName, samples_use = NULL,
                       fun = st_any_pred, return_df = FALSE, ...) {
    cg <- colGeometry(x,
        type = colGeometryName, sample_id = samples_use,
        withDimnames = TRUE
    )
    cg$barcode <- rownames(cg)
    cg <- cg[, c("barcode", "geometry")]
    cgs <- split(cg, colData(x)$sample_id[colData(x)$sample_id %in% samples_use])
    out <- lapply(names(cgs), function(s) {
        if ("sample_id" %in% names(y)) {
            y_use <- y[y$sample_id == s, ]
        } else {
            y_use <- y
        }
        o <- fun(cgs[[s]], y_use, ...)
        if (!return_df) {
            names(o) <- cgs[[s]]$barcode
        } else {
            rownames(o) <- cgs[[s]]$barcode
        }
        o
    })
    if (!return_df) {
        out <- Reduce(c, out)
        out <- out[rownames(cg)]
    } else {
        out <- do.call(rbind, out)
        out <- out[rownames(cg), ]
        out$barcode <- NULL
    }
    out
}

#' Binary predicates for geometry of each cell/spot and annotation
#'
#' This function finds binary predicates for the geometry of each cell/spot
#' (i.e. \code{colGeometry}) and an annotation geometry for each sample. For
#' example, whether each Visium spot intersects with the tissue boundary in each
#' sample.
#'
#' @param sfe An SFE object.
#' @param colGeometryName Name of column geometry for the predicate.
#' @param annotGeometryName Name of annotation geometry for the predicate.
#' @param sample_id Which sample(s) to operate on. Can be "all" to indicate all
#'   samples.
#' @param pred Predicate function to use, defaults to
#'   \code{\link{st_intersects}}.
#' @return For \code{annotPred}, a logical vector of the same length as the
#'   number of columns in the sample(s) of interest, with barcodes (or
#'   corresponding column names of sfe) as names. For \code{annotNPred}, a
#'   numeric vector of the same length as the number of columns in the sample(s)
#'   of interest with barcodes as names, indicating the number of geometries
#'   in the \code{annotGeometry} of interest returns TRUE for the predicate for
#'   each each geometry in the \code{colGeometry} of interest.
#' @concept Geometric operations
#' @export
#' @seealso annotOp
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' # Whether each spot is in tissue
#' in_tissue <- annotPred(sfe, "spotPoly", annotGeometryName = "tissueBoundary")
#' # How many nuclei are there in each Visium spot
#' n_nuclei <- annotNPred(sfe, "spotPoly", annotGeometryName = "nuclei")
annotPred <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                      sample_id = "all", pred = st_intersects) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
    .annot_fun(sfe, ag,
        colGeometryName = colGeometryName,
        samples_use = sample_id, fun = st_any_pred, pred = pred
    )
}

#' @rdname annotPred
#' @export
annotNPred <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                       sample_id = "all", pred = st_intersects) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
    .annot_fun(sfe, ag,
        colGeometryName = colGeometryName,
        samples_use = sample_id, fun = st_n_pred, pred = pred
    )
}

#' Binary operations for geometry of each cell/spot and annotation
#'
#' Just like \code{\link{annotPred}}, but performs the operation rather than
#' predicate. For example, this function would return the geometry of the
#' intersections between each Visium spot and the tissue boundary for each
#' sample, rather than whether each Visium spot intersects the tissue boundary.
#' In case one cell/spot gets broken up into multiple geometries, the union of
#' those geometries will be taken, so each cell/spot will only get one geometry.
#'
#' @inheritParams annotPred
#' @param op A binary operation function for the geometries. Defaults to
#' \code{\link{st_intersection}}.
#' @return A \code{sf} data frame with \code{geometry} column containing the
#' geometries and corresponding column names of sfe as row names. There is no
#' guarantee that the returned geometries are valid or preserve the geometry
#' class (e.g. when the intersection of polygons result into a line of a point).
#' @concept Geometric operations
#' @export
#' @seealso annotPred
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' # Get the intersection of myofibers with each Visium spot
#' myofibers_on_spots <- annotOp(sfe, "spotPoly",
#'     annotGeometryName = "myofiber_simplified"
#' )
annotOp <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                    sample_id = "all", op = st_intersection) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
    cg <- colGeometry(sfe, type = colGeometryName, sample_id = sample_id)
    sample_index <- colData(sfe)$sample_id %in% sample_id
    out <- .crop_geometry(cg, ag,
        samples_use = sample_id, op = op,
        sample_col = colData(sfe)$sample_id[sample_index],
        id_col = "barcode"
    )
    if (!"sample_id" %in% names(cg)) out$sample_id <- NULL
    out <- out[rownames(cg), ]
    rownames(out) <- rownames(cg)
    out
}

.annot_summ <- function(x, y, pred, summary_fun) {
    out <- st_drop_geometry(st_join(x, y, join = pred))
    barcodes <- out$barcode
    out <- out[, names(out) != "barcode", drop = FALSE]
    aggregate(out, list(barcode = barcodes),
        FUN = summary_fun, simplify = TRUE,
        drop = FALSE
    )
}

#' Summarize attributes of an annotGeometry for each cell/spot
#'
#' In SFE objects, the annotation geometries don't have to correspond to the
#' dimensions of the gene count matrix, so there generally is no one to one
#' mapping between annotation geometries and cells/spots. However, it may be
#' interesting to relate attributes of annotation geometries to cell/spots so
#' the attributes can be related to gene expression. This function summarizes
#' attributes of an \code{annotGeometry} for each cell/spot by a geometric
#' predicate with a \code{colGeometry}.
#'
#' @inheritParams annotPred
#' @param annotColNames Character, column names of the \code{annotGeometry} of
#'   interest, to indicate the columns to summarize. Columns that are absent
#'   from the \code{annotGeometry} are removed. The column cannot be "geometry"
#'   or "barcode".
#' @param summary_fun Function for the summary, defaults to \code{mean}.
#' @return A data frame whose row names are the relevant column names of
#'   \code{sfe}, and each column of which is the summary of each column
#'   specified in \code{annotColName}.
#' @concept Geometric operations
#' @importFrom sf st_join
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' s <- annotSummary(sfe, "spotPoly", "myofiber_simplified",
#'     annotColNames = c("area", "convexity")
#' )
annotSummary <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                         annotColNames = 1L, sample_id = "all",
                         pred = st_intersects, summary_fun = mean) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
    cols_use <- intersect(
        annotColNames,
        names(ag)[!names(ag) %in% c("barcode", "geometry")]
    )
    if (!length(cols_use)) {
        stop(
            "None of the columns specified in annotColNames is present in the",
            " annotGeometry ", annotGeometryName
        )
    }
    ag <- ag[, cols_use] # Geometry is sticky
    .annot_fun(sfe, ag,
        colGeometryName = colGeometryName,
        samples_use = sample_id, fun = .annot_summ, pred = pred,
        summary_fun = summary_fun, return_df = TRUE
    )
}

.valid_bbox <- function(bbox) {
    # Only for one vector
    if (bbox["xmin"] > bbox["xmax"] || bbox["ymin"] > bbox["ymax"])
        stop("Min limit is greater than max limit")
}

.check_bbox <- function(bbox) {
    if (!is.numeric(bbox))
        stop("bbox must be a numeric vector or matrix.")
    if (is.vector(bbox)) {
        if (!setequal(names(bbox), c("xmin", "xmax", "ymin", "ymax")))
            stop("Vector bbox must be a vector of length 4 with names xmin, xmax, ymin, ymax in any order.")
        .valid_bbox(bbox)
    } else if (is.matrix(bbox)) {
        if (!setequal(rownames(bbox), c("xmin", "xmax", "ymin", "ymax"))) {
            stop("Matrix bbox must have rownames xmin, xmax, ymin, ymax, in any order.")
        }
        if (is.null(colnames(bbox))) {
            stop("Matrix bbox must have colnames to indicate sample ID.")
        }
        apply(bbox, 2, .valid_bbox)
    }
}

.bbox2sf <- function(bbox, sample_id) {
    if (is.vector(bbox) || is(bbox, "bbox")) {
        bbox <- matrix(bbox, ncol = 1, dimnames = list(names(bbox), sample_id[1]))
    } else {
        samples_use <- intersect(colnames(bbox), sample_id)
        if (!length(samples_use))
            stop("No bounding boxes for samples specified.")
        bbox <- bbox[, samples_use, drop = FALSE]
    }
    bbox <- bbox[c("xmin", "ymin", "xmax", "ymax"),, drop = FALSE]
    bbox_sfc <- apply(bbox, 2, function(x) st_as_sfc(st_bbox(x)),
                      simplify = FALSE)
    bbox_sfc <- Reduce(c, bbox_sfc)
    bbox_sf <- st_sf(geometry = bbox_sfc)
    if (!is.null(sample_id)) bbox_sf$sample_id <- colnames(bbox)
    bbox_sf
}

#' Crop an SFE object with a geometry
#'
#' Returns an SFE object whose specified \code{colGeometry} returns \code{TRUE}
#' with a geometric predicate function (usually intersects) with another
#' geometry of interest. This can be used to subset an SFE object with a tissue
#' boundary or histological region polygon, or crop away empty spaces. After
#' cropping, not only will the cells/spots be subsetted, but also all geometries
#' will be cropped.
#'
#' 3D geometries are allowed, but geometric operations can only be performed in
#' x and y but not z.
#'
#' @param x An SFE object.
#' @param y An object of class \code{sf}, \code{sfg}, \code{sfc} with which to
#'   crop the SFE object, or a bounding box with the format of the output of
#'   \code{\link{bbox,SpatialFeatureExperiment-method}}.
#' @param colGeometryName Column geometry to used to indicate which cells/spots
#'   to keep.
#' @param sample_id Samples to crop. Optional when only one sample is present.
#'   Can be multiple samples, or "all", which means all samples. For multiple
#'   samples, \code{sf} data frame \code{y} may have column \code{sample_id}
#'   indicating which geometry subsets which sample or matrix \code{y} may
#'   indicate sample specific bounding boxes in its column names. Only samples
#'   included in the indicated sample IDs are subsetted. If sample is not
#'   indicated in \code{y}, then the same geometry or bounding box is used to
#'   subset all samples specified in the \code{sample_id} argument.
#' @param op A geometric operation function to crop the geometries in the SFE
#'   object. Only \code{\link{st_intersection}} and \code{\link{st_difference}}
#'   are allowed. If "intersection", then only things inside \code{y} is kept
#'   after cropping. If "difference", then only things outside \code{y} is kept.
#' @param keep_whole Character vector, can be one or more of "col" and "annot"
#'   to keep whole items from \code{colGeometries} or \code{annotGeometries},
#'   keeping geometries that partially intersect with \code{y} whole. This can
#'   greatly speed up code while not breaking geometries into multiple pieces.
#'   Can also be "none" so all geometries are actually cropped.
#' @param cover Logical, whether the geometries in \code{x} must be entirely
#'   covered by \code{y} if \code{op = st_intersection} or whether \code{x} must
#'   be entirely outside \code{y} if \code{op = st_difference}. Only relevant
#'   when \code{keep_whole != "none"}.
#' @return An SFE object. There is no guarantee that the geometries after
#'   cropping are still all valid or preserve the original geometry class.
#' @concept Geometric operations
#' @importFrom sf st_intersection st_union st_agr st_covered_by st_disjoint
#' @importFrom lifecycle deprecated is_present deprecate_warn
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' # Subset sfe to only keep spots on tissue
#' sfe_on_tissue <- crop(sfe, tissueBoundary(sfe),
#'     colGeometryName = "spotPoly",
#'     sample_id = "Vis5A"
#' )
crop <- function(x, y = NULL, colGeometryName = 1L, sample_id = "all",
                 op = st_intersection,
                 keep_whole = "none",
                 cover = FALSE) {
    if (!(identical(op, sf::st_intersection) || identical(op, sf::st_difference))) {
        stop("op must be either st_intersection or st_difference.")
    }
    keep_whole <- match.arg(keep_whole, choices = c("none", "col", "annot"),
                            several.ok = TRUE)
    if (length(keep_whole) > 1L && "none" %in% keep_whole) {
        keep_whole <- keep_whole[keep_whole != "none"]
    }
    is_difference <- identical(op, sf::st_difference)
    pred <- if (is_difference) {
        if (cover) st_disjoint else function(x, y, sparse = TRUE) !st_covered_by(x, y, sparse = sparse)
    } else {
        if (cover) st_covered_by else st_intersects
    }
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (!is(y, "sf") && !is(y, "sfc") && !is(y, "sfg")) {
        # y should be bbox, either named vector or matrix with samples in columns
        .check_bbox(y)
        if (is.matrix(y) && is.null(sample_id)) sample_id <- colnames(y)
        y <- .bbox2sf(y, sample_id)
    }
    if ("sample_id" %in% names(y)) {
        samples_use <- intersect(unique(y$sample_id), sample_id)
        if (!length(samples_use))
            stop("No bounding boxes for samples specified.")
    } else {
        samples_use <- sample_id
    }
    preds <- .annot_fun(x, y, colGeometryName,
        samples_use = samples_use,
        fun = st_any_pred, pred = pred
    )
    # Don't remove anything from other samples
    other_bcs <- setdiff(colnames(x), names(preds))
    other_samples <- setNames(rep(TRUE, length(other_bcs)), other_bcs)
    preds <- c(preds, other_samples)
    preds <- preds[colnames(x)]
    out <- x[, preds]
    samples_rmed <- setdiff(sample_id, sampleIDs(out))
    if (!any(preds)) {
        warning("No cells/spots left after cropping")
    } else if (length(samples_rmed)) {
        warning(
            "Sample(s) ", paste(samples_rmed, collapse = ", "),
            " was/were removed by the cropping."
        )
    }
    if (!"col" %in% keep_whole) {
        # colGeometries should have already been subsetted here
        # Also actually crop all geometries for the samples of interest.
        colGeometries(out) <- lapply(colGeometries(out), .crop_geometry,
                                     y = y,
                                     samples_use = samples_use, op = op,
                                     id_col = "barcode",
                                     sample_col = colData(out)$sample_id
        )
    }
    if ("annot" %in% keep_whole) {
        annotGeometries(out) <- lapply(annotGeometries(out), function(x) {
            inds <- x$sample_id %in% sample_id
            other_samples <- x[!inds,]
            in_sample <- x[inds,]
            preds <- pred(in_sample, y, sparse = FALSE)
            if (!any(preds)) other_samples else
                rbind(in_sample[preds,], other_samples)
        })
    } else {
        annotGeometries(out) <- lapply(annotGeometries(out), .crop_geometry,
                                       y = y,
                                       samples_use = samples_use, op = op,
                                       remove_empty = TRUE
        )
    }
    if (length(rowGeometries(out))) {
        for (s in samples_use) {
            rgs <- rowGeometries(out, sample_id = s)
            rowGeometries(out, sample_id = s) <-
                lapply(rgs, function(r) {
                    r <- .crop_geometry(r, y = y, op = op, sample_col = s,
                                        samples_use = s)
                })
        }
    }
    out <- .crop_imgs(out, bbox(out, sample_id = samples_use, include_image = FALSE))
    out
}

.bbox_sample_g <- function(gs) {
    # Assume that it's already properly subsetted for the sample
    # For one sample_id, multiple geometries
    bboxes <- lapply(gs, st_bbox)
    aggBboxes(bboxes)
}

.bbox_sample <- function(sfe, sample_id, include_images = FALSE,
                         include_row = FALSE) { # For one sample
    sample_index <- colData(sfe)$sample_id == sample_id
    cgs <- as.list(int_colData(sfe)[["colGeometries"]][sample_index, ,
        drop = FALSE
    ])
    cgs_bboxes <- lapply(cgs, st_bbox)
    ags <- annotGeometries(sfe)
    if (length(ags)) {
        ags <- lapply(ags, function(g) g[g$sample_id == sample_id, ])
        ags_bboxes <- lapply(ags, st_bbox)
    } else {
        ags_bboxes <- NULL
    }

    if (include_row) {
        rgs <- rowGeometries(sfe)
        if (length(rgs) && length(sampleIDs(sfe)) > 1L) {
            pattern <- paste0(sample_id, "$")
            has_sample <- grepl(pattern, names(rgs))
            rgs <- rgs[has_sample]
        }
        if (length(rgs)) rgs_bboxes <- lapply(rgs, st_bbox) else rgs_bboxes <- NULL
    } else rgs_bboxes <- NULL
    if (isEmpty(imgData(sfe))) include_images <- FALSE
    if (include_images) {
        img_exts <- lapply(imgData(sfe)$data, ext)
    } else img_exts <- NULL
    all_bboxes <- c(cgs_bboxes, ags_bboxes, rgs_bboxes, img_exts)
    aggBboxes(all_bboxes)
}

#' Find bounding box of SFE objects
#'
#' Find bounding box of the union of all \code{colGeometries} and
#' \code{annotGeometries} of each sample in the SFE object. This can be used to
#' remove empty space so the tissue and geometries have one corner at the origin
#' so all samples will be on comparable coordinates.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @param sample_id Sample(s) whose bounding box(es) to find. The bounding box
#'   would be for the union of all \code{colGeometries} and
#'   \code{annotGeometries} associated with each sample.
#' @param include_images Logical, whether the bounding boxes should include
#'   image extents. Defaults to \code{FALSE} because often the image has a lot
#'   of empty space surrounding the tissue.
#' @param include_row Logical, whether the bounding boxes should include
#' \code{rowGeometries}, defaults to \code{TRUE}.
#' @return For one sample, then a named vector with names \code{xmin},
#'   \code{ymin}, \code{xmax}, and \code{ymax} specifying the bounding box. For
#'   multiple samples, then a matrix whose columns are samples and whose rows
#'   delineate the bounding box.
#' @concept Geometric operations
#' @aliases bbox
#' @importFrom sf st_bbox
#' @importFrom S4Vectors isEmpty
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' bbox(sfe, sample_id = "Vis5A")
setMethod("bbox", "SpatialFeatureExperiment", function(sfe, sample_id = "all",
                                                       include_images = FALSE,
                                                       include_row = TRUE) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    out <- lapply(sample_id, .bbox_sample, sfe = sfe,
                  include_images = include_images,
                  include_row = include_row)
    # Just as in sf for empty geometries
    if (!length(out))
        return(c(xmin = NA, ymin = NA, xmax = NA, ymax = NA))
    if (length(out) == 1L) {
        out <- out[[1]]
    } else {
        out <- do.call(cbind, out)
        colnames(out) <- sample_id
    }
    out
})

#' Affine transfortaion of SFE object in histological space
#'
#' These functions perform affine transformations on SFE objects, including all
#' geometries and images. The transformation is performed on each sample as a
#' whole. This differs from functions such as \code{\link{mirrorImg}} in that
#' \code{mirrorImg} and \code{rotateImg} transform the image with the center at
#' the center of the image itself. In contrast, the center of transformation
#' here is the center of the bounding box of the entire sample, including
#' images.
#'
#' For images that are not loaded into memory, \code{rotateImg} will load
#' \code{SpatRasterImage} into memory and all image operations except translate
#' will load \code{BioFormatsImage} into memory.
#'
#' @inheritParams terra::flip
#' @inheritParams transposeImg
#' @inheritParams scaleImg
#' @inheritParams affineImg
#' @inheritParams rotateImg
#' @param sfe An SFE object.
#' @param sample_id Sample(s) to transform.
#' @param maxcell Rotating \code{SpatRasterImage} will convert it into
#'   \code{ExtImage}, loading the image into memory. This argument specifies the
#'   maximum number of pixels in the image loaded into memory. The image will be
#'   down sampled to approximately this number of pixels.
#' @param v Vector to spatially translate the SFE object.
#' @return An SFE object with the sample(s) transformed.
#' @name SFE-transform
#' @concept Geometric operations
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' sfe2 <- transpose(sfe)
#' sfe3 <- mirror(sfe)
#'
NULL

#' @rdname SFE-transform
#' @export
transpose <- function(sfe, sample_id = "all", maxcell = NULL, filename = "") {
    .transform_samples(sfe, sample_id, geometry_fun = .transpose_geometry,
                       img_fun = .transpose_img, maxcell = maxcell,
                       filename = filename)
}

#' @rdname SFE-transform
#' @export
mirror <- function(sfe, sample_id = "all",
                   direction = c("vertical", "horizontal"), maxcell = NULL,
                   filename = "") {
    .transform_samples(sfe, sample_id, geometry_fun = .mirror_geometry,
                       img_fun = .mirror_img, maxcell = maxcell,
                       filename = filename, direction = direction)
}

#' @rdname SFE-transform
#' @export
rotate <- function(sfe, sample_id = "all", degrees, maxcell = 1e7) {
    .transform_samples(sfe, sample_id, geometry_fun = .rotate_geometry,
                       img_fun = .rotate_img, maxcell = maxcell, degrees = degrees)
}

#' @rdname SFE-transform
#' @export
translate <- function(sfe, sample_id = "all", v) {
    .transform_samples(sfe, sample_id, geometry_fun = .translate_geometry,
                       img_fun = translateImg, v = v, use_bbox = FALSE)
}

#' @rdname SFE-transform
#' @export
scale <- function(sfe, sample_id = "all", factor) {
    .transform_samples(sfe, sample_id, geometry_fun = .scale_geometry,
                       img_fun = .scale_ext, factor = factor)
}

#' @rdname SFE-transform
#' @export
affine <- function(sfe, sample_id = "all", M, v, maxcell = 1e7) {
    .transform_samples(sfe, sample_id, geometry_fun = .affine_geometry,
                       img_fun = affineImg, M = M, v = v, maxcell = maxcell,
                       use_bbox = FALSE)
}

#' Remove empty space
#'
#' For each sample independently, all geometries and \code{spatialCoords} are
#' translated so the origin is at the minimum coordinates of the bounding box
#' of all geometries of the sample. This way coordinates of different samples
#' will be more comparable. This removes empty space in the images if present.
#'
#' @param sfe An SFE object.
#' @param sample_id Sample to remove empty space.
#' @return An SFE object with empty space removed.
#' @note Unlike other functions in this package, this function operates on all
#' samples by default.
#' @concept Geometric operations
#' @export
#' @examples
#' library(SFEData)
#' library(SingleCellExperiment)
#' sfe <- McKellarMuscleData("full")
#' # Only keep spots on tissue
#' sfe <- sfe[, colData(sfe)$in_tissue]
#' # Move the coordinates of the tissue
#' sfe <- removeEmptySpace(sfe)
removeEmptySpace <- function(sfe, sample_id = "all") {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    bboxes <- bbox(sfe, sample_id, include_images = FALSE)
    sfe <- .crop_imgs(sfe, bboxes)
    if (length(sample_id) == 1L) {
        sfe <- translate(sfe, sample_id, v = -bboxes[c("xmin", "ymin")])
    } else {
        for (s in sample_id) {
            sfe <- translate(sfe, s, v = -bboxes[c("xmin", "ymin"), s])
        }
    }
    if (!is.null(int_metadata(sfe)$orig_bbox)) {
        int_metadata(sfe)$orig_bbox <- bboxes
    } else {
        if (!is.matrix(bboxes)) {
            bboxes <- matrix(bboxes, ncol = 1)
            rownames(bboxes) <- c("xmin", "ymin", "xmax", "ymax")
            colnames(bboxes) <- sample_id
        }
        int_metadata(sfe)$orig_bbox <- cbind(int_metadata(sfe)$orig_bbox, bboxes)
    }
    sfe
}
