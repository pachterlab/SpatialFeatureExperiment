#' Images in SpatialFeatureExperiment object
#'
#' \code{SpatialFeatureExperiment} and the \code{Voyager} package work with
#' images differently from \code{SpatialExperiment}. In SFE and
#' \code{Voyager}'s, plotting functions for SFE objects, the images are read
#' with \code{\link{terra::rast}} and represented as \code{SpatRaster}, so the
#' image is not entirely loaded into memory unless necessary. Plotting will not
#' load a large image into memory; rather the image will be downsampled and the
#' downsampled version is plotted.
#'
#' @param extent A numeric vector of length 4 with names of the set xmin, ymin,
#'   xmax, and ymax, specifying the extent of the image.
#' @param scale_fct Scale factor -- multiply pixel coordinates in full
#'   resolution image by this scale factor should yield pixel coordinates in a
#'   different resolution. \code{extent} takes precedence over \code{scale_fct}.
#' @note If the image is already a GeoTIFF file that already has an extent, then
#'   the extent associated with the file will be honored and the \code{extent}
#'   and \code{scale_fct} arguments are ignored.
#' @return Methods for \code{SpatRasterImage} return a modified
#'   \code{SpatRasterImage}, and methods for SFE return a modified SFE object.
#' @importClassesFrom SpatialExperiment VirtualSpatialImage
#' @importFrom SpatialExperiment addImg mirrorImg imgData imgData<- imgRaster
#'   imgSource
#' @importFrom terra ext ext<-
#' @importClassesFrom terra SpatRaster
#' @name SFE-image
NULL

#' @rdname SFE-image
#' @export
setClass("SpatRasterImage", contains="VirtualSpatialImage",
         slots=c(image="SpatRaster"))

#' @rdname SFE-image
#' @export
setMethod("addImg", "SpatialFeatureExperiment",
          function(x, file, sample_id, image_id, extent = NULL,
                   scale_fct = 1) {
              file <- normalizePath(file, mustWork = TRUE)
              sample_id <- .check_sample_id(x, sample_id)
              if (!is.null(extent)) {
                  if (!is.numeric(extent))
                      stop("extent must be a numeric vector.")
                  if (!setequal(names(extent), c("xmin", "xmax", "ymin", "ymax")))
                      stop("Names of extent must be of the set 'xmin', 'xmax', 'ymin', and 'ymax'.")
              } else {
                  if (!(scale_fct > 0 && is.numeric(scale_fct))) {
                      stop("scale_fct must be a positive number.")
                  }
              }
              # From SPE
              # check that image entry doesn't already exist
              idx <- tryCatch(
                  error=function(e) e,
                  SpatialExperiment:::.get_img_idx(x, sample_id, image_id))

              if (!inherits(idx, "error"))
                  stop("'imgData' already contains an entry with",
                       sprintf(
                           " 'image_id = %s' and 'sample_id = %s'",
                           dQuote(image_id), dQuote(sample_id)))

              # get & add valid 'imgData' entry
              df <- .get_imgData(file, sample_id, image_id, extent, scale_fct)
              imgData(x) <- rbind(imgData(x), df)
              return(x)
          })

.get_imgData <- function(file, sample_id, image_id, extent = NULL,
                         scale_fct = 1, flip = FALSE) {
    # What if extent is already present?
    w <- tryCatch(rast(file), warning = function(w) w)
    if (is(w, "warning")) {
        # No extent in tif file
        suppressWarnings(img <- rast(file))
        if (!is.null(extent)) ext(img) <- extent[c("xmin", "xmax", "ymin", "ymax")]
        else {
            ext(img) <- as.vector(ext(img)) / scale_fct
        }
    } else {
        # Honor extent in the GeoTIFF file
        img <- w
        flip <- FALSE
    }
    if (flip) img <- terra::flip(img)
    spi <- new("SpatRasterImage", image = img)
    DataFrame(
        sample_id,
        image_id,
        data=I(list(spi)),
        scaleFactor=scale_fct)
}

# Just to make the SPE functions happy
#' @rdname SFE-image
#' @export
setMethod("transposeImg", "SpatRasterImage",
          function(x, left = TRUE) {
              x@image <- terra::trans(x@image)
              x
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "SpatRasterImage",
          function(x, direction = "vertical") {
              x@image <- terra::flip(x@image, direction = direction)
              x
          })

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "SpatialFeatureExperiment",
          function(x, sample_id = NULL, image_id = NULL) {
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  new <- lapply(old, transposeImg)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "SpatialFeatureExperiment",
          function(x, sample_id=NULL, image_id=NULL, direction = "vertical") {
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  new <- lapply(old, mirrorImg, direction = direction)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

#' @rdname SFE-image
#' @export
setMethod("imgRaster", "SpatRasterImage", function(x) x@image)

#' @rdname SFE-image
#' @export
setMethod("imgSource",
          "SpatRasterImage",
          function(x, path=FALSE) {
              NA_character_
          })

.crop_imgs <- function(x, bboxes) {
    # Only works for SpatRaster
    if (nrow(imgData(x))) {
        samples <- sampleIDs(x)
        if (length(samples) == 1L) {
            bboxes <- matrix(bboxes, ncol = 1)
            colnames(bboxes) <- samples
        }
        new_imgs <- lapply(samples, function(s) {
            img_data <- imgData(x)$data[imgData(x)$sample_id == s]
            bbox_use <- ext(bboxes[c("xmin", "xmax", "ymin", "ymax"),s])
            lapply(img_data, function(img) {
                new("SpatRasterImage", image = terra::crop(img@image, bbox_use))
            })
        })
        new_imgs <- unlist(new_imgs, recursive = FALSE)
        imgData(x)$data <- I(new_imgs)
    }
    x
}
