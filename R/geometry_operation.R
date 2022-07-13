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
#' @param pred A geometric binary predicate function, such as \code{\link{st_intersects}}.
#' It should return an object of class \code{sgbp}, for sparse predicates.
#' @return For \code{st_any_*}, a logical vector indicating whether each geometry in \code{x} intersects
#' (or other predicates such as is covered by) anything in \code{y}. Simplified
#' from the \code{sgbp} results which indicate which item in \code{y} each item
#' in \code{x} intersects, which might not always be relevant. For \code{st_n_*},
#' an integer vector indicating the number of geometries in y returns TRUE for
#' each geometry in x.
#' @export
#' @importFrom sf st_intersects st_agr<- st_drop_geometry st_as_sfc st_cast
#' st_is_empty st_disjoint
#' @importFrom stats aggregate
#' @examples
#' library(sf)
#' pts <- st_sfc(st_point(c(.5,.5)), st_point(c(1.5, 1.5)), st_point(c(2.5, 2.5)))
#' pol <- st_polygon(list(rbind(c(0,0), c(2,0), c(2,2), c(0,2), c(0,0))))
#' st_any_pred(pts, pol, pred = st_disjoint)
#' st_any_intersects(pts, pol)
#' st_n_pred(pts, pol, pred = st_disjoint)
#' st_n_intersects(pts, pol)
st_any_pred <- function(x, y, pred) lengths(pred(x, y)) > 0L

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
  # g has column sample_id as in colGeometry and annotGeometry
  # g should also have row names if it is colGeometry
  if (!id_col %in% names(g)) {
    rm_id <- TRUE
    if (!is.null(rownames(g))) {
      g[[id_col]] <- rownames(g)
    } else {
      g[[id_col]] <- as.character(seq_len(nrow(g)))
    }
  } else rm_id <- FALSE
  if (!is.null(sample_col)) g$sample_id <- sample_col
  gs <- split(g, g$sample_id)
  gs_sub <- lapply(names(gs), function(s) {
    if (s %in% samples_use) {
      if ("sample_id" %in% names(y))
        y_use <- st_geometry(y[y$sample_id == s,])
      else y_use <- st_geometry(y)
      .g <- gs[[s]][,c("geometry", id_col)]
      st_agr(.g) <- "constant"
      o <- op(.g, y_use)
      # Aggregate in case cropping broke some items into multiple pieces
      if (any(!rownames(o) %in% rownames(.g))) {
        o <- aggregate(o, by = setNames(list(id = o[[id_col]]), id_col),
                       FUN = unique)
      }
      merge(o, st_drop_geometry(gs[[s]]), by = id_col, all = TRUE)
    }
    else gs[[s]]
  })
  gs_sub <- do.call(rbind, gs_sub)
  # Convert st_GEOMETRY to a more specific type
  if (st_geometry_type(gs_sub, by_geometry = FALSE) == "GEOMETRY") {
    gs_sub <- st_cast(gs_sub)
  }
  gs_sub <- gs_sub[,names(g)]
  rownames(gs_sub) <- rownames(g)[match(gs_sub[[id_col]], g[[id_col]])]
  if (rm_id) gs_sub[[id_col]] <- NULL
  if (remove_empty) gs_sub <- gs_sub[!st_is_empty(gs_sub),]
  gs_sub
}

.annot_fun <- function(x, y, colGeometryName, samples_use = NULL, fun = st_any_pred,
                       return_df = FALSE, ...) {
  cg <- colGeometry(x, type = colGeometryName, sample_id = samples_use,
                    withDimnames = TRUE)
  cg$barcode <- rownames(cg)
  cg <- cg[, c("barcode", "geometry")]
  cgs <- split(cg, colData(x)$sample_id[colData(x)$sample_id %in% samples_use])
  out <- lapply(names(cgs), function(s) {
    if ("sample_id" %in% names(y))
      y_use <- y[y$sample_id == s,]
    else y_use <- y
    o <- fun(cgs[[s]], y_use, ...)
    if (!return_df) names(o) <- cgs[[s]]$barcode
    else rownames(o) <- cgs[[s]]$barcode
    o
  })
  if (!return_df) {
    out <- Reduce(c, out)
    out <- out[rownames(cg)]
  } else {
    out <- do.call(rbind, out)
    out <- out[rownames(cg),]
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
                      sample_id = NULL, pred = st_intersects) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
  .annot_fun(sfe, ag, colGeometryName = colGeometryName, samples_use = sample_id,
             fun = st_any_pred, pred = pred)
}

#' @rdname annotPred
#' @export
annotNPred <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                      sample_id = NULL, pred = st_intersects) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
  .annot_fun(sfe, ag, colGeometryName = colGeometryName, samples_use = sample_id,
             fun = st_n_pred, pred = pred)
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
#' @export
#' @seealso annotPred
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' # Get the intersection of myofibers with each Visium spot
#' myofibers_on_spots <- annotOp(sfe, "spotPoly", annotGeometryName = "myofiber_simplified")
annotOp <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                    sample_id = NULL, op = st_intersection) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
  cg <- colGeometry(sfe, type = colGeometryName, sample_id = sample_id)
  out <- .crop_geometry(cg, ag, samples_use = sample_id, op = op,
                        sample_col = colData(sfe)$sample_id[colData(sfe)$sample_id %in% sample_id],
                        id_col = "barcode")
  if (!"sample_id" %in% names(cg)) out$sample_id <- NULL
  out <- out[rownames(cg),]
  rownames(out) <- rownames(cg)
  out
}

.annot_summ <- function(x, y, pred, summary_fun) {
  out <- st_drop_geometry(st_join(x, y, join = pred))
  barcodes <- out$barcode
  out <- out[,names(out) != "barcode", drop = FALSE]
  aggregate(out, list(barcode = barcodes), FUN = summary_fun, simplify = TRUE, drop = FALSE)
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
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' s <- annotSummary(sfe, "spotPoly", "myofiber_simplified",
#'                   annotColNames = c("area", "convexity"))
annotSummary <- function(sfe, colGeometryName = 1L, annotGeometryName = 1L,
                         annotColNames = 1L, sample_id = NULL,
                         pred = st_intersects, summary_fun = mean) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  ag <- annotGeometry(sfe, type = annotGeometryName, sample_id = sample_id)
  cols_use <- intersect(annotColNames, names(ag)[!names(ag) %in% c("barcode", "geometry")])
  if (!length(cols_use))
    stop("None of the columns specified in annotColNames is present in the annotGeometry ",
         annotGeometryName)
  ag <- ag[,cols_use] # Geometry is sticky
  .annot_fun(sfe, ag, colGeometryName = colGeometryName, samples_use = sample_id,
             fun = .annot_summ, pred = pred, summary_fun = summary_fun,
             return_df = TRUE)
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
#' @param x An SFE object.
#' @param y An object of class \code{sf}, \code{sfg}, or \code{sfc} with which
#'   to crop the SFE object. Optional if \code{xmin}, \code{xmax}, \code{ymin},
#'   and \code{ymax} are specified for a bounding box.
#' @param colGeometryName Column geometry to used to indicate which cells/spots
#'   to keep.
#' @param sample_id Samples to crop. Optional when only one sample is present.
#'   Can be multiple samples, or "all", which means all samples. For multiple
#'   samples, \code{y} may have column \code{sample_id} indicating which
#'   geometry subsets which sample. Only samples included in the
#'   \code{sample_id} column are subsetted. If there is no \code{sample_id}
#'   column or \code{y} is not specified, then the same geometry or bounding box
#'   is used to subset all samples specified in the \code{sample_id} argument.
#' @param pred A geometric binary predicate function to indicate which
#'   cells/spots to keep, defaults to \code{\link{st_intersects}}.
#' @param op A geometric operation function to crop the geometries in the SFE
#'   object. Defaults to \code{\link{st_intersection}}.
#' @param xmin Minimum x coordinate of bounding box. Ignored if \code{y} is
#'   specified.
#' @param xmax Maximum x coordinate of bounding box.
#' @param ymin Minimum y coordinate of bounding box.
#' @param ymax Maximum y coordinate of bounding box.
#' @return An SFE object. There is no guarantee that the geometries after
#'   cropping are still all valid or preserve the original geometry class.
#' @importFrom sf st_intersection st_union st_agr
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' # Subset sfe to only keep spots on tissue
#' sfe_on_tissue <- crop(sfe, tissueBoundary(sfe), colGeometryName = "spotPoly",
#'                       sample_id = "Vis5A")
#' # Subset sfe to only keep what's within a bounding box
#' # All geometries will be cropped
#' # sample_id is optional when only one sample is present
#' sfe_cropped <- crop(sfe, colGeometryName = "spotPoly",
#'                     xmin = 5500, xmax = 6500, ymin = 13500, ymax = 14500)
crop <- function(x, y = NULL, colGeometryName = 1L, sample_id = NULL,
                 pred = st_intersects, op = st_intersection, xmin = NULL,
                 xmax = NULL, ymin = NULL, ymax = NULL) {
  sample_id <- .check_sample_id(x, sample_id, one = FALSE)
  if (is.null(y)) {
    y <- st_sf(geometry = st_as_sfc(st_bbox(c(xmin = xmin, xmax = xmax,
                                              ymin = ymin, ymax = ymax),
                                            crs = NA)),
               sf_column_name = "geometry")
  }
  if (!is(y, "sf") && !is(y, "sfc") && !is(y, "sfg")) {
    stop("y must be a sf, sfc, or sfg object")
  }
  if ("sample_id" %in% names(y)) {
    samples_use <- intersect(unique(y$sample_id), sample_id)
  } else {
    samples_use <- sample_id
  }
  preds <- .annot_fun(x, y, colGeometryName, samples_use = samples_use, fun = st_any_pred,
                      pred = pred)
  # Don't remove anything from other samples
  other_bcs <- setdiff(colnames(x), names(preds))
  other_samples <- setNames(rep(TRUE, length(other_bcs)), other_bcs)
  preds <- c(preds, other_samples)
  preds <- preds[colnames(x)]
  out <- x[, preds]
  # colGeometries should have already been subsetted here
  # Also actually crop all geometries for the samples of interest. Leave rowGeometry alone for now.
  colGeometries(out) <- lapply(colGeometries(out), .crop_geometry, y = y,
                               samples_use = samples_use, op = op,
                               id_col = "barcode", sample_col = colData(out)$sample_id)
  annotGeometries(out) <- lapply(annotGeometries(out), .crop_geometry, y = y,
                                 samples_use = samples_use, op = op,
                                 remove_empty = TRUE)
  samples_rmed <- setdiff(sample_id, sampleIDs(out))
  if (length(samples_rmed)) {
    warning("Sample(s) ", paste(samples_rmed, collapse = ", "),
            " was/were removed by the cropping.")
  }
  out
}

.bbox_sample <- function(sfe, sample_id) {
  cgs <- as.list(int_colData(sfe)[["colGeometries"]][colData(sfe)$sample_id == sample_id,,drop = FALSE])
  cgs_bboxes <- lapply(cgs, st_bbox)
  ags <- annotGeometries(sfe)
  if (length(ags)) {
    ags <- lapply(ags, function(g) g[g$sample_id == sample_id,])
    ags_bboxes <- lapply(ags, st_bbox)
    all_bboxes <- c(cgs_bboxes, ags_bboxes)
  } else all_bboxes <- cgs_bboxes
  all_bboxes <- do.call(rbind, all_bboxes)
  c(xmin = min(all_bboxes[,"xmin"], na.rm = TRUE),
    ymin = min(all_bboxes[,"ymin"], na.rm = TRUE),
    xmax = max(all_bboxes[,"xmax"], na.rm = TRUE),
    ymax = max(all_bboxes[,"ymax"], na.rm = TRUE))
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
#' @return For one sample, then a named vector with names \code{xmin},
#'   \code{ymin}, \code{xmax}, and \code{ymax} specifying the bounding box. For
#'   multiple samples, then a matrix whose columns are samples and whose rows
#'   delineate the bounding box.
#' @aliases bbox
#' @importFrom sf st_bbox
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' bbox(sfe, sample_id = "Vis5A")
setMethod("bbox", "SpatialFeatureExperiment", function(sfe, sample_id = NULL) {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  out <- lapply(sample_id, .bbox_sample, sfe = sfe)
  if (length(out) == 1L) {
    out <- out[[1]]
  } else {
    out <- do.call(cbind, out)
    colnames(out) <- sample_id
  }
  out
})

.translate <- function(sfe, sample_id, v) {
  # spatialCoords
  inds <- colData(sfe)$sample_id == sample_id
  spatialCoords(sfe)[inds,] <- spatialCoords(sfe)[inds,] + v
  # colGeometries
  for (n in colGeometryNames(sfe)) {
    cg <- colGeometry(sfe, n, sample_id = sample_id)
    cg$geometry <- cg$geometry + v
    colGeometry(sfe, n, sample_id = sample_id) <- cg
  }
  # annotGeometries
  if (!is.null(annotGeometries(sfe))) {
    for (n in annotGeometryNames(sfe)) {
      ag <- annotGeometry(sfe, n, sample_id = sample_id)
      ag$geometry <- ag$geometry + v
      annotGeometry(sfe, n, sample_id = sample_id) <- ag
    }
  }
  sfe
}

# Burning question: What to do with rowGeometry? I suppose I can require that
# the geometries are specified separately for each sample. I'll deal with that
# later. Perhaps much later, as I don't see it relevant to Visium.

#' Remove empty space
#'
#' For each sample independently, all geometries and \code{spatialCoords} are
#' translated so the origin is at the minimum coordinates of the bounding box
#' of all geometries of the sample. This way coordinates of different samples
#' will be more comparable.
#'
#' @param sfe An SFE object.
#' @param sample_id Sample to remove empty space.
#' @return An SFE object with empty space removed.
#' @note Unlike other functions in this package, this function operates on all
#' samples by default.
#' @export
#' @examples
#' library(SFEData)
#' library(SingleCellExperiment)
#' sfe <- McKellarMuscleData("full")
#' # Only keep spots on tissue
#' sfe <- sfe[,colData(sfe)$in_tissue]
#' # Move the coordinates of the tissue
#' sfe <- removeEmptySpace(sfe)
removeEmptySpace <- function(sfe, sample_id = "all") {
  sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
  bboxes <- bbox(sfe, sample_id)
  if (length(sample_id) == 1L) {
    sfe <- .translate(sfe, sample_id, v = -bboxes[c("xmin", "ymin")])
  } else {
    for (s in sample_id) {
      sfe <- .translate(sfe, s, v = -bboxes[c("xmin", "ymin"), s])
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
