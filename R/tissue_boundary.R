.get_mask_boundary <- function(mask, n_pieces = 1, simplify = TRUE, dTolerance = 0) {
    mask <- toSpatRasterImage(mask, save_geotiff = FALSE)
    names(mask) <- "value"
    polys <- terra::as.polygons(mask)
    polys <- polys[polys$value==1,]
    polys <- st_as_sf(polys) |> st_geometry()
    polys <- st_cast(polys, "POLYGON")
    area <- st_area(polys)
    if (n_pieces < length(area)) 
        polys <- polys[order(area, decreasing = TRUE)[seq_len(n_pieces)]]
    if (simplify) {
        #polys <- rmapshaper::ms_simplify(polys, ...)
        polys <- st_simplify(polys, dTolerance = dTolerance)
    }
    polys
}

.get_factor <- function(d) {
    # for clahe
    f <- as.integer(gmp::factorize(d))
    f <- f[f <= 15]
    cp <- cumprod(f)
    max(cp[cp <= 15])
}
#' @importFrom EBImage normalize otsu closing makeBrush fillHull
.get_tb1 <- function(sfe, sample_id, image_id, 
                     image_type, channel, n_pieces, resolution, maxcell,
                     fill_holes, simplify, dTolerance) {
    # For one sample
    img <- getImg(sfe, image_id = image_id, sample_id = sample_id)
    if (inherits(img, "BioFormatsImage")) 
        img <- toExtImage(img, resolution = resolution, channel = channel)
    else if (inherits(img, "SpatRasterImage"))
        img <- toExtImage(img, maxcell = maxcell, channel = channel)
    # rowSums is a bit confusing; with dims = 2 it sums all channels, i.e. the 3rd dimension
    if (length(dim(img)) > 2L) img <- rowSums(img, dims = 2)
    img <- EBImage::normalize(img)
    if (image_type == "brightfield") {
        mask <- img < otsu(img)
    } else {
        fctx <- .get_factor(dim(img)[1])
        fcty <- .get_factor(dim(img)[2])
        img <- clahe(img, nx = fctx, ny = fcty)
        mask <- img > otsu(img)
        mask <- closing(mask, makeBrush(7))
    }
    if (fill_holes) mask <- fillHull(mask)
    .get_mask_boundary(mask, n_pieces = n_pieces, simplify = simplify,
                       dTolerance = dTolerance)
}

#' Get tissue boundary from histology image
#'
#' This function gets the tissue boundary from image and makes sure that it is
#' properly aligned with the geometries in the SFE object.
#'
#' @param sfe An SFE object with images
#' @param sample_id Sample id(s) whose tissue boundaries are to be found.
#' @param image_type Character, either "brightfield" or "fluorescent"
#' @param channel Channel to use for tissue segmentation. If \code{NULL} use
#'   average of all channels.
#' @param resolution Integer, which resolution to use for tissue boundary in a
#'   pyramidal OME-TIFF stack. Only applicable to \code{BioFormatsImage}. Note
#'   that the image will be loaded into memory and you usually don't need the
#'   highest resolution for the tissue boundary.
#' @param maxcell Max number of pixels when loading \code{SpatRasterImage} into
#'   memory.
#' @param n_pieces Number of pieces; only this number of largest pieces are
#'   kept. Smaller pieces will be considered debris and removed. Can be a vector
#'   of the same length as \code{sample_id}, or if \code{sample_id} is "all"
#'   then same length as the number of samples in the SFE object to specify
#'   different number of pieces in different samples. If \code{n_pieces} is
#'   length 1 while there are multiple samples, then the same number is applied
#'   to all samples.
#' @param fill_holes Logical, whether to fill holes in the tissue, to only get
#'   the outer outline.
#' @param simplify Logical, whether to simplify the output polygon.
#' @param dTolerance Distance tolerance when simplifying the polygon, in the
#'   same unit as the geometries in the SFE object.
#' @return A \code{sf} data frame with columns \code{sample_id} and
#'   \code{geometry}.
#' @importFrom sf st_simplify
#' @export
#' 
getTissueBoundaryImg <- function(sfe, sample_id = NULL, image_id = NULL,
                              image_type = c("brightfield", "fluorescent"),
                              channel = NULL, n_pieces = 1, resolution = 4,
                              maxcell = 1e7,
                              fill_holes = FALSE, simplify = TRUE, 
                              dTolerance = 0) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    image_id <- image_id %||% imageIDs(sfe)[1]
    if (length(n_pieces) != 1 && length(n_pieces) != length(sample_id)) {
        stop("n_pieces must be of either length 1 or same length as sample_id.")
    }
    image_type <- match.arg(image_type)
    if (image_type == "fluorescent") check_installed("gmp")
    if (length(n_pieces) > 1L) {
        og <- mapply(.get_tb1, sample_id = sample_id, n_pieces = n_pieces,
                     MoreArgs = list(sfe = sfe, image_type = image_type,
                                     channel = channel, resolution = resolution,
                                     maxcell = maxcell, fill_holes = fill_holes,
                                     simplify = simplify, dTolerance = dTolerance),
                     SIMPLIFY = FALSE)
        og <- st_sfc(unlist(og, recursive = FALSE))
    } else {
        og <- .get_tb1(sfe, sample_id, image_id, image_type, channel, n_pieces,
                       resolution, maxcell, fill_holes, simplify, dTolerance)
    }
    st_sf(sample_id = sample_id, geometry = og, agr = "constant", crs = NA)
}

.get_tb1_concave <- function(sfe, sample_id, colGeometryName, ratio, allow_holes) {
    cs <- st_union(colGeometry(sfe, colGeometryName, sample_id = sample_id))
    st_concave_hull(cs, ratio = ratio, allow_holes = allow_holes)
}

#' Get tissue boundary from concave hull of cell geometries
#'
#' The concave hull will be a smoothed outline, and the smoothness can be
#' adjusted with the \code{ratio} parameter. Run \code{\link{findDebrisCells}}
#' to remove small bits outside the main piece tissue before running this
#' function.
#' 
#' @inheritParams getTissueBoundaryImg
#' @inheritParams sf::st_concave_hull
#' @return A \code{sf} data frame with columns \code{sample_id} and
#'   \code{geometry}.
#' @export
getTissueBoundaryConcave <- function(sfe, sample_id = NULL, colGeometryName = 1L,
                                     ratio = 0.01, allow_holes = TRUE) {
    sample_id <- .check_sample_id(sfe, sample_id, one = FALSE)
    og <- lapply(sample_id, .get_tb1_concave, sfe = sfe, 
                 colGeometryName = colGeometryName, ratio = ratio, 
                 allow_holes = allow_holes)
    og <- st_sfc(unlist(og, recursive = FALSE))
    st_sf(sample_id = sample_id, geometry = og, agr = "constant", crs = NA)
}
