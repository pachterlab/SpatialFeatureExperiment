#' Images in SpatialFeatureExperiment object
#'
#' \code{SpatialFeatureExperiment} and the \code{Voyager} package work with
#' images differently from \code{SpatialExperiment}. In SFE and
#' \code{Voyager}'s, plotting functions for SFE objects, the images are read
#' with \code{\link{rast}} and represented as \code{SpatRaster}, so the image is
#' not entirely loaded into memory unless necessary. Plotting will not load a
#' large image into memory; rather the image will be downsampled and the
#' downsampled version is plotted.
#'
#' @inheritParams terra::flip
#' @param sample_id Which sample the image is associated with. Use
#'   \code{\link{sampleIDs}} to get sample IDs present in the SFE object.
#' @param image_id Image ID, such as "lowres" and "hires" for Visium data and
#' "DAPI" and "PolyT" for Vizgen MERFISH data.
#' @param file File from which to read the image.
#' @param extent A numeric vector of length 4 with names of the set xmin, ymin,
#'   xmax, and ymax, specifying the extent of the image.
#' @param scale_fct Scale factor -- multiply pixel coordinates in full
#'   resolution image by this scale factor should yield pixel coordinates in a
#'   different resolution. \code{extent} takes precedence over \code{scale_fct}.
#' @note If the image is already a GeoTIFF file that already has an extent, then
#'   the extent associated with the file will be honored and the \code{extent}
#'   and \code{scale_fct} arguments are ignored. Also, when the image is
#'   transposed, it is flipped about the axis going from top left to bottom
#'   right.
#' @return Methods for \code{SpatRasterImage} return a modified
#'   \code{SpatRasterImage}, and methods for SFE return a modified SFE object.
#' @importClassesFrom SpatialExperiment VirtualSpatialImage
#' @importFrom SpatialExperiment addImg mirrorImg imgData imgData<- imgRaster
#'   imgSource getImg
#' @importFrom terra ext ext<-
#' @importClassesFrom terra SpatRaster
#' @importClassesFrom EBImage Image
#' @importFrom methods setClassUnion
#' @concept Image and raster
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' img_path <- system.file(file.path("extdata", "sample01", "outs", "spatial",
#' "tissue_lowres_image.png"),
#' package = "SpatialFeatureExperiment")
#' sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
#' scale_fct = 0.023)
#' img <- getImg(sfe)
#' # SpatRasterImage method
#' img_t <- transposeImg(img)
#' # SFE method
#' sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "lowres")
#' @name SFE-image
#' @aliases transposeImg
NULL

# Not exported by terra
setClass("PackedSpatRaster",
         representation (
             definition = "character",
             values = "matrix",
             attributes = "list"
         ),
         prototype (
             attributes = list()
         )
)

setClassUnion("GeneralizedSpatRaster", members = c("SpatRaster", "PackedSpatRaster"))

#' @rdname SFE-image
#' @export
setClass("SpatRasterImage", contains="VirtualSpatialImage",
         slots=c(image="GeneralizedSpatRaster"))

#' @rdname SFE-image
#' @export
setClass("BioFormatsImage", contains = "VirtualSpatialImage",
         slots = c(path = "character", ext = "numeric"))
#' @rdname SFE-image
#' @export
setClass("EBImage", contains = "VirtualSpatialImage",
         slots = c(image = "Image", ext = "numeric"))

#' @rdname SFE-image
#' @export
SpatRasterImage <- function(img) {
    new("SpatRasterImage", image = img)
}

#' @rdname SFE-image
#' @export
EBImage <- function(img, ext) {
    new("EBImage", image = img, ext = ext)
}

#' @rdname SFE-image
#' @export
BioFormatsImage <- function(path, ext) {
    new("BioFormatsImage", path = path, ext = ext)
}

#' @rdname SFE-image
#' @export
setMethod("ext", "BioFormatsImage", function(x) x@ext)

#' @rdname SFE-image
#' @export
setMethod("ext", "EBImage", function(x) x@ext)

#' @rdname SFE-image
#' @export
setMethod("ext", "SpatRasterImage", function(x) ext(x@image))

.set_ext <- function(x, value) {
    x@ext <- value
    x
}

#' @rdname SFE-image
#' @export
setReplaceMethod("ext", c("BioFormatsImage", "numeric"), .set_ext)

#' @rdname SFE-image
#' @export
setReplaceMethod("ext", c("EBImage", "numeric"), .set_ext)

#' @rdname SFE-image
#' @export
setReplaceMethod("ext", c("SpatRasterImage", "numeric"),
                 function(x, value) {
                     ext(x@image) <- value
                     x
                 })

setValidity("BioFormatsImage", function(x) {
    outs <- character(2)
    if (!file.exists(x@path)) outs[1] <- paste0("File ", x@path, " does not exist.")
    out[2] <- tryCatch(.check_bbox(x@ext), error = function(e) e$message)
    out <- out[out != ""]
    if (length(out)) return(out) else TRUE
})

setValidity("EBImage", function(x) {
    out <- tryCatch(.check_bbox(x@ext), error = function(e) e$message)
    if (length(out)) return(out) else TRUE
})

# How useful is scaleFactor here?
toEBImage <- function(bfi, resolution = 4L, scaleFactor = NULL) {
    check_installed(c("xml2", "RBioFormats"))
    # PhysicalSizeX seems to be a standard field
    if (length(resolution) != 1L ||
        !isTRUE(all.equal(floor(resolution), resolution))) {
        stop("resolution must be integer of length 1.")
    }
    if (!is.null(scaleFactor) &&
        (length(scaleFactor) != 1L || !is.numeric(scaleFactor) ||
         any(scaleFactor <= 0L))) {
        stop("scaleFactor must be a positive number.")
    }
    file <- imgSource(bfi)
    xml_meta <- RBioFormats::read.omexml(file) |>
        xml2::read_xml() |> xml2::as_list()
    psx <- attr(xml_meta$OME$Image$Pixels$TiffData, "PhysicalSizeX")
    psy <- attr(xml_meta$OME$Image$Pixels$TiffData, "PhysicalSizeY")
    if (is.null(psx) && is.null(psy) && is.null(scaleFactor)) {
        warning("Physical pixel size absent from image metadata. ",
                "Unspecified scaleFactor assumed to be 1.")
    }
    if (!is.null(psx) && is.null(psy)) {
        message("Assuming equal physical pixel size on x and y.")
        psy <- psx
    } else if (is.null(psx) && !is.null(psy)) {
        message("Assuming equal physical pixel size on x and y.")
        psx <- psy
    }
    if (!is.null(psx) && !is.null(psy) && !is.null(scaleFactor)) {
        sfx <- 1/psx
        sfy <- 1/psy # I think x and y sizes are usually the same but just in case
        if (!isTRUE(all.equal(sfx, scaleFactor)))
            message("Physical pixel size found in metadata. Ignoring scaleFactor.")
    }
    if (is.null(psx) && is.null(psy) && !is.null(scaleFactor)) {
        sfx <- sfy <- scaleFactor
    }
    if (!is.na(ext(bfi))) {
        bbox_use <- ext(bfi)
        x_nms <- c("xmin", "xmax")
        y_nms <- c("ymin", "ymax")
        # Convert to full res pixel space
        bbox_use[x_nms] <- bbox_use[x_nms] * sfx
        bbox_use[y_nms] <- bbox_use[y_nms] * sfy
        # Convert to lower res pixel space
        if (resolution > 1L) {
            meta <- RBioFormats::read.metadata(file) |>
                RBioFormats::coreMetadata()
            sizeX <- meta[[resolution]]$sizeX
            sizeY <- meta[[resolution]]$sizeY
            sizeX_full <- meta[[1]]$sizeX
            sizeY_full <- meta[[1]]$sizeY
            sfx2 <- sizeX/sizeX_full
            sfy2 <- sizeY/sizeY_full
            bbox_use[x_nms] <- bbox_use[x_nms] * sfx2
            bbox_use[y_nms] <- bbox_use[y_nms] * sfy2
        }
        min_nms <- c("xmin", "ymin")
        max_nms <- c("xmax", "ymax")
        bbox_use[min_nms] <- floor(bbox_use[min_nms])
        bbox_use[max_nms] <- ceiling(bbox_use[max_nms])
        subset_use <- list(x = seq(bbox_use["xmin"], bbox_use["xmax"], by = 1L),
                           y = seq(bbox_use["ymin"], bbox_use["ymax"], by = 1L))
    } else subset_use <- list()

    img <- RBioFormats::read.image(file = file,
                                   resolution = resolution,
                                   filter.metadata = TRUE,
                                   read.metadata = FALSE,
                                   normalize = FALSE,
                                   subset = subset_use)
    EBImage(img, ext(bfi))
}

toSpatRasterImage <- function(bfi, resolution = 4L, scaleFactor = NULL) {
    check_installed("RBioFormats")
    # Only for OME-TIFF, haven't tested on other BioFormats
    img <- toEBImage(bfi, resolution, scaleFactor)
    img_fn <- gsub("ome.", "", imgSource(bfi))
    message(">>> Saving lower resolution images with `.tif` (non OME-TIFF) format:",
            paste0("\n", img_fn))
    RBioFormats::write.image(img@image, file = img_fn, force = TRUE)
    SpatRasterImage(rast(img_fn))
}

#' @rdname SFE-image
#' @export
setMethod("addImg", "SpatialFeatureExperiment",
          function(x, imageSource, sample_id, image_id,
                   extent = NULL, scale_fct = 1, file = deprecated()) {
              if (lifecycle::is_present(file)) {
                  deprecate_warn("1.6.0", "SpatialFeatureExperiment::addImg(file = )",
                                 "SpatialFeatureExperiment::addImg(imageSource = )")
                  imageSource <- file
              }
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
              df <- .get_imgData(imageSource, sample_id, image_id, extent, scale_fct)
              imgData(x) <- rbind(imgData(x), df)
              return(x)
          })

.get_imgData <- function(img, sample_id, image_id, extent = NULL,
                         scale_fct = 1, flip = FALSE) {
    # EBImage
    if (is(img, "Image")) {
        spi <- EBImage(img)
    } else if (.path_valid2(img)) {
        e <- tryCatch(rast(img), error = function(e) e)
        if (is(e, "error")) {
            spi <- BioFormatsImage(img, extent)
            if (flip) spi <- mirrorImg(spi, direction = "vertical")
        } else {
            # What if extent is already present?
            w <- tryCatch(rast(img), warning = function(w) w)
            if (is(w, "warning")) {
                # No extent in tif file
                suppressWarnings(img <- rast(img))
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
        }
    } else {
        spi <- SpatialExperiment:::.get_imgData(img, scale_fct, sample_id, image_id)
    }

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
          function(x, ...) {
              x@image <- terra::trans(imgRaster(x))
              x
          })

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "BioFormatsImage",
          function(x, resolution = 4L, scaleFactor = NULL) {
              x <- toEBImage(x, resolution, scaleFactor)
              transposeImg(x)
          })

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "EBImage",
          function(x, ...) {
              x@image <- EBImage::transpose(x@image)
              x
          })

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "SpatialFeatureExperiment",
          function(x, sample_id = NULL, image_id = NULL,
                   resolution = 4L) {
              sample_id <- .check_sample_id(x, sample_id, one = TRUE)
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  sfct <- imgData(x)$scaleFactor[idx]
                  new <- mapply(transposeImg, x = old, scaleFactor = sfct,
                                MoreArgs = list(resolution = resolution))
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "SpatRasterImage",
          function(x, direction = c("vertical", "horizontal"), ...) {
              direction <- match.arg(direction)
              x@image <- terra::flip(imgRaster(x), direction = direction)
              x
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "BioFormatsImage",
          function(x, direction = c("vertical", "horizontal"),
                   resolution = 4L, scaleFactor = NULL) {
              direction <- match.arg(direction)
              x <- toEBImage(x, resolution, scaleFactor)
              mirrorImg(x, direction)
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "EBImage",
          function(x, direction = c("vertical", "horizontal"), ...) {
              direction <- match.arg(direction)
              fun <- if (direction == "vertical") EBImage::flip else EBImage::flop
              x@image <- fun(x@image)
              x
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "SpatialFeatureExperiment",
          function(x, sample_id=NULL, image_id=NULL, direction = "vertical",
                   resolution = 4L) {
              sample_id <- .check_sample_id(x, sample_id, one = TRUE)
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  sfct <- imgData(x)$scaleFactor[idx]
                  new <- mapply(mirrorImg, x = old, scaleFactor = sfct,
                                MoreArgs = list(resolution = resolution,
                                                direction = direction))
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })
# Just realized that SpatialExperiment only rotates by multiples of 90 degrees
# Can apply that and reconstruct rast. Will load into memory.
#' @rdname SFE-image
#' @export
setMethod("rotateImg", "SpatRasterImage", # Deal with rotating SpatRaster?
          function(x, degree, ...) {
              # Not sure what exactly to do. I think convert to EBImage as well.
              # TODO: write function to convert SPI to EBI
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "BioFormatsImage",
          function(x, degree, resolution = 4L, scaleFactor = NULL) {
              x <- toEBImage(x, resolution, scaleFactor)
              rotateImg(x, degree)
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "EBImage",
          function(x, degree, ...) {
              x@image <- EBImage::rotate(x@image, degree)
              x
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "SpatialFeatureExperiment",
          function(x, sample_id=NULL, image_id=NULL, degree,
                   resolution = 4L) {
              sample_id <- .check_sample_id(x, sample_id, one = TRUE)
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  sfct <- imgData(x)$scaleFactor[idx]
                  new <- mapply(rotateImg, x = old, scaleFactor = sfct,
                                MoreArgs = list(resolution = resolution,
                                                degree = degree))
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

#' @rdname SFE-image
#' @export
setMethod("imgRaster", "SpatRasterImage", function(x) {
    if (is(x@image, "PackedSpatRaster")) unwrap(x@image)
    else x@image
})

#' @rdname SFE-image
#' @export
setMethod("imgRaster", "BioFormatsImage", function(x, resolution = 4L) {
    toEBImage(x, resolution) |> imgRaster()
})

#' @rdname SFE-image
#' @export
setMethod("imgRaster", "EBImage", function(x) x@image)

#' @rdname SFE-image
#' @export
#' @importFrom terra sources
setMethod("imgSource",
          "SpatRasterImage",
          function(x) {
              sources(imgRaster(x))
          })

#' @rdname SFE-image
#' @export
setMethod("imgSource", "BioFormatsImage", function(x) x@path)

# TODO: S4 crop methods for BioFormatsImage, EBImage, and LoadedSpatialImage
# EBImage crops with matrix-like subsetting of pixels
setMethod("crop", "BioFormatsImage", function(x, bbox) {

})
.crop_imgs <- function(x, bboxes) {
    # Only works for SpatRaster
    if (nrow(imgData(x))) {
        samples <- sort(sampleIDs(x))
        imgData(x) <- imgData(x)[order(imgData(x)$sample_id),]
        if (length(samples) == 1L) {
            bboxes <- matrix(bboxes, ncol = 1, dimnames = list(names(bboxes), samples))
        }
        new_imgs <- lapply(samples, function(s) {
            img_data <- imgData(x)$data[imgData(x)$sample_id == s]
            bbox_use <- ext(bboxes[c("xmin", "xmax", "ymin", "ymax"),s])
            lapply(img_data, function(img) {
                new("SpatRasterImage", image = terra::crop(imgRaster(img), bbox_use, snap = "out"))
                # Deal with BioFormatsImage
            })
        })
        new_imgs <- unlist(new_imgs, recursive = FALSE)
        imgData(x)$data <- I(new_imgs)
    }
    x
}
