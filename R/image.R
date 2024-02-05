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
#'   imgSource getImg rotateImg
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
setMethod("ext", "SpatRasterImage", function(x) ext(x@image) |> as.vector())

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
                     ext(x@image) <- value[c("xmin", "xmax", "ymin", "ymax")]
                     x
                 })

setValidity("BioFormatsImage", function(object) {
    outs <- character(2)
    if (!file.exists(object@path)) outs[1] <- paste0("File ", object@path, " does not exist.")
    out[2] <- tryCatch(.check_bbox(object@ext), error = function(e) e$message)
    out <- out[out != ""]
    if (length(out)) return(out) else TRUE
})

setValidity("EBImage", function(object) {
    out <- tryCatch(.check_bbox(object@ext), error = function(e) e$message)
    if (length(out)) return(out) else TRUE
})

.toEBImage <- function(x, resolution = 4L) {
    check_installed(c("xml2", "RBioFormats"))
    # PhysicalSizeX seems to be a standard field
    if (length(resolution) != 1L ||
        !isTRUE(all.equal(floor(resolution), resolution))) {
        stop("resolution must be integer of length 1.")
    }
    file <- imgSource(x)
    xml_meta <- RBioFormats::read.omexml(file) |>
        xml2::read_xml() |> xml2::as_list()
    psx <- attr(xml_meta$OME$Image$Pixels$TiffData, "PhysicalSizeX")
    psy <- attr(xml_meta$OME$Image$Pixels$TiffData, "PhysicalSizeY")
    if (is.null(psx) && is.null(psy)) {
        stop("Physical pixel size absent from image metadata.")
    }
    if (!is.null(psx) && is.null(psy)) {
        message("Assuming equal physical pixel size on x and y.")
        psy <- psx
    } else if (is.null(psx) && !is.null(psy)) {
        message("Assuming equal physical pixel size on x and y.")
        psx <- psy
    }
    sfx <- 1/psx
    sfy <- 1/psy
    meta <- RBioFormats::read.metadata(file) |>
        RBioFormats::coreMetadata()
    # Only 1 resolution is present
    one_res <- "sizeX" %in% names(meta)
    if (one_res) {
        resolution <- 1L
        sizeX_full <- meta$sizeX
        sizeY_full <- meta$sizeY
    } else {
        sizeX_full <- meta[[1]]$sizeX
        sizeY_full <- meta[[1]]$sizeY
    }
    min_nms <- c("xmin", "ymin")
    max_nms <- c("xmax", "ymax")
    x_nms <- c("xmin", "xmax")
    y_nms <- c("ymin", "ymax")
    if (!is.na(ext(x))) {
        bbox_use <- ext(x)
        # Convert to full res pixel space
        bbox_use[x_nms] <- bbox_use[x_nms] * sfx
        bbox_use[y_nms] <- bbox_use[y_nms] * sfy
        # Convert to lower res pixel space
        if (resolution > 1L) {
            sizeX <- meta[[resolution]]$sizeX
            sizeY <- meta[[resolution]]$sizeY
            sfx2 <- sizeX/sizeX_full
            sfy2 <- sizeY/sizeY_full
            bbox_use[x_nms] <- bbox_use[x_nms] * sfx2
            bbox_use[y_nms] <- bbox_use[y_nms] * sfy2
        } else sfx2 <- sfy2 <- 1
        bbox_use[min_nms] <- floor(bbox_use[min_nms])
        bbox_use[max_nms] <- ceiling(bbox_use[max_nms])
        subset_use <- list(x = seq(bbox_use["xmin"], bbox_use["xmax"], by = 1L),
                           y = seq(bbox_use["ymin"], bbox_use["ymax"], by = 1L))
        # Extent should account for the pixels
        ext_use <- bbox_use
        ext_use[x_nms] <- ext_use[x_nms]/(sfx*sfx2)
        ext_use[y_nms] <- ext_use[y_nms]/(sfy*sfy2)
    } else {
        ext_use <- c(xmin = 0, ymin = 0, xmax = sizeX_full/sfx, ymax = sizeY_full/sfy)
        subset_use <- list()
    }

    img <- RBioFormats::read.image(file = file,
                                   resolution = resolution,
                                   filter.metadata = TRUE,
                                   read.metadata = FALSE,
                                   normalize = FALSE,
                                   subset = subset_use)
    EBImage(img, ext_use)
}

# SpatRasterImage method: allow downsampling first
# From tidyterra and ggspavis::plotVisium
.resample_spat <- function(r, maxcell = 50000) {
    if (terra::ncell(r) > 1.1 * maxcell) {
        r <- terra::spatSample(r, maxcell,
                               as.raster = TRUE,
                               method = "regular"
        )
    }
    return(r)
}

# TODO: refactor functions related to geom_spi_rgb in Voyager since some of the code got moved here
.toEBImage2 <- function(x, maxcell = 1e7) {
    # 1e7 comes from the number of pixels in resolution = 4L in the ome.tiff
    if (length(names(x)) == 3L) {
        names(x) <- c("r", "g", "b")
        # Remove RGB settings, better plot without it
        terra::RGB(x) <- NULL
    }
    x <- .resample_spat(x, maxcell)
    # Problem: in EBImage, the 1st dimension is x/column in the image
    if (length(names(x)) == 3L)
        out <- terra::as.array(x) |> aperm(c(2,1,3)) |> Image(colormode = "Color")
    else
        out <- terra::as.array(x)[,,1] |> t() |> Image(colormode = "Grayscale")
    out
}

#' Convert images to EBImage
#'
#' The \code{EBImage} class is a thin wrapper around the \code{Image} class in
#' \code{EBImage} so it inherits from \code{VirtualSpatialImage} as required by
#' \code{SpatialExperiment} and has extent as used in Voyager's plotting
#' functions. This function converts \code{SpatRasterImage} (thin wrapper around
#' the class in \code{terra}) and \code{BioFormatsImage} into \code{EBImage} for
#' image operations as implemented in the \code{EBImage} package.
#'
#' @param x Either a \code{BioFormatsImage} or \code{SpatRasterImage} object.
#' @param resolution Integer, which resolution in the \code{BioFormatsImage} to
#' read and convert. Defaults to 4, which is a lower resolution.
#' @return A \code{EBImage} object. The image is loaded into memory.
#' @name toEBImage
#' @aliases toEBImage
#' @export
setMethod("toEBImage", "BioFormatsImage", .toEBImage)

#' @rdname toEBImage
#' @export
setMethod("toEBImage", "SpatRasterImage", .toEBImage2)

toSpatRasterImage <- function(x, resolution = 4L) {
    check_installed("RBioFormats")
    # Only for OME-TIFF, haven't tested on other BioFormats
    img <- toEBImage(x, resolution)
    img_fn <- gsub("ome.", "", imgSource(x))
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

.get_imgData <- function(img, sample_id, image_id, extent = NA,
                         scale_fct = 1, flip = FALSE) {
    # EBImage
    if (is(img, "Image")) {
        spi <- EBImage(img, extent)
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
                if (!is.na(extent)) ext(img) <- extent[c("xmin", "xmax", "ymin", "ymax")]
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
              # What terra does to extent: swap xmin and xmax with ymin and ymax
              x
          })

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "BioFormatsImage",
          function(x, resolution = 4L) {
              x <- toEBImage(x, resolution)
              transposeImg(x)
          })

.trans_extent <- function(e) {
    ext_use <- e[c("xmin", "xmax", "ymin", "ymax")]
    names(ext_use) <- c("ymin", "ymax", "xmin", "xmax")
    ext_use
}

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "EBImage",
          function(x, ...) {
              x@image <- EBImage::transpose(x@image)
              # Extent
              ext(x) <- .trans_extent(ext(x))
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
                  new <- lapply(old, transposeImg, resolution = resolution)
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
                   resolution = 4L) {
              direction <- match.arg(direction)
              x <- toEBImage(x, resolution)
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
                  new <- lapply(old, mirrorImg, resolution = resolution,
                                direction = direction)
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })
# Just realized that SpatialExperiment only rotates by multiples of 90 degrees
# Can apply that and reconstruct rast. Will load into memory.
#' @rdname SFE-image
#' @export
setMethod("rotateImg", "SpatRasterImage", # Deal with rotating SpatRaster?
          function(x, degree, maxcell = 1e7, ...) {
              # Not sure what exactly to do. I think convert to EBImage as well.
              stopifnot(
                  length(degrees) == 1,
                  is.numeric(degrees),
                  degrees %% 90 == 0)
              x <- toEBImage(x, maxcell)
              rotateImg(x, degree)
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "BioFormatsImage",
          function(x, degree, resolution = 4L, ...) {
              # Only allow multiples of 90 degree, deal with extent
              stopifnot(
                  length(degrees) == 1,
                  is.numeric(degrees),
                  degrees %% 90 == 0)
              x <- toEBImage(x, resolution)
              rotateImg(x, degree)
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "EBImage",
          function(x, degree, ...) {
              stopifnot(
                  length(degrees) == 1,
                  is.numeric(degrees),
                  degrees %% 90 == 0)
              x@image <- EBImage::rotate(x@image, degree)
              # Deal with extent
              if (degrees %% 180 > 0) {
                  ext(x) <- .trans_extent(ext(x))
              }
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
                  new <- lapply(old, rotateImg, resolution = resolution,
                                maxcell = maxcell,
                                degree = degree)
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

# TODO: imgRaster setter function since here I want to allow image processing
# like adjusting brightness and contrast, blurring, sharpening, opening, closing, and so on
# but what if it changes the extent?

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


setMethod("cropImg", "SpatRasterImage", function(x, bbox) {
    x@image <- terra::crop(x@image, bbox, snap = "out")
    x
})

setMethod("cropImg", "BioFormatsImage", function(x, bbox) {
    bbox_old <- ext(x)
    bbox_use <- c(xmin = max(bbox_old["xmin"], bbox["xmin"]),
                  xmax = min(bbox_old["xmax"], bbox["xmax"]),
                  ymin = max(bbox_old["ymin"], bbox["ymin"]),
                  ymax = min(bbox_old["ymax"], bbox["ymax"]))
    ext(x) <- bbox_use
    x
})

setMethod("cropImg", "EBImage", function(x, bbox) {
    # Convert bbox to pixel range based on ext(x)
    bbox_old <- ext(x)
    dim_old <- dim(x@image)
    sfx <- dim_old[1]/(bbox_old["xmax"] - bbox_old["xmin"])
    sfy <- dim_old[2]/(bbox_old["ymax"] - bbox_old["ymin"])
    bbox_use <- bbox
    bbox_use[c("xmin", "xmax")] <- bbox_old[c("xmin", "xmax")] * sfx
    bbox_use[c("ymin", "ymax")] <- bbox_old[c("ymin", "ymax")] * sfy
    bbox_use[c("xmin", "ymin")] <- floor(bbox_use[c("xmin", "ymin")])
    bbox_use[c("xmax", "ymax")] <- ceiling(bbox_use[c("xmax", "ymax")])

    if (length(dim_old) == 3L)
        x@image <- x@image[seq(bbox_use["xmin"], bbox_use["xmax"]),
                           seq(bbox_use["ymin"], bbox_use["ymax"]),]
    else
        x@image <- x@image[seq(bbox_use["xmin"], bbox_use["xmax"]),
                           seq(bbox_use["ymin"], bbox_use["ymax"])]

    bbox_new <- bbox_use
    bbox_new[c("xmin", "xmax")] <- bbox_new[c("xmin", "xmax")] / sfx
    bbox_new[c("ymin", "ymax")] <- bbox_new[c("ymin", "ymax")] / sfy
    ext(x) <- bbox_new
    x
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
