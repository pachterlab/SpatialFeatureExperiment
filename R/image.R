# SpatRasterImage====================

#' SpatRaster representation of images in SFE objects
#'
#' \code{SpatialFeatureExperiment} and the \code{Voyager} package work with
#' images differently from \code{SpatialExperiment}. In SFE and
#' \code{Voyager}'s, plotting functions for SFE objects, the images can be read
#' with \code{\link{rast}} and represented as \code{SpatRaster}, so the image is
#' not entirely loaded into memory unless necessary. Plotting will not load a
#' large image into memory; rather the image will be downsampled and the
#' downsampled version is plotted.
#'
#' @param img A \code{\link{SpatRaster}} or \code{\link{PackedSpatRaster}} object.
#' @note If the image is already a GeoTIFF file that already has an extent, then
#'   the extent associated with the file will be honored and the \code{extent}
#'   and \code{scale_fct} arguments are ignored. Also, when the image is
#'   transposed, it is flipped about the axis going from top left to bottom
#'   right.
#' @return Methods for \code{SpatRasterImage} return a modified
#'   \code{SpatRasterImage}, and methods for SFE return a modified SFE object.
#' @importClassesFrom SpatialExperiment VirtualSpatialImage
#' @importFrom SpatialExperiment addImg mirrorImg imgData imgData<- imgRaster
#'   imgSource getImg rotateImg rmvImg
#' @importFrom terra ext ext<-
#' @importClassesFrom terra SpatRaster
#' @importClassesFrom EBImage Image
#' @importFrom methods setClassUnion
#' @concept Image and raster
#' @examples
#' # Example code
#' @name SpatRasterImage
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

#' @rdname SpatRasterImage
#' @export
setClass("SpatRasterImage", contains="VirtualSpatialImage",
         slots=c(image="GeneralizedSpatRaster"))

#' @rdname SpatRasterImage
#' @export
SpatRasterImage <- function(img) {
    new("SpatRasterImage", image = img)
}

# BioFormatsImage====================

#' On disk representation of BioFormats images in SFE object
#'
#' At present, the \code{BioFormatsImage} is designed for OME-TIFF from Xenium
#' and has not been tested on other formats that can be read with
#' \code{BioFormats}. The image is not loaded into memory, and when it is, the
#' the \code{BioFormatsImage} object is converted into \code{\link{EBImage}}
#' because the loaded image is of a class that inherits from
#' \code{\link{EBImage::Image}}. The \code{\link{EBImage}} class is a thin
#' wrapper inheriting from \code{VirtualSpatialImage} so it's compatible with
#' \code{SpatialExperiment} from which SFE is derived.
#'
#' @param path Path to an OME-TIFF image file.
#' @param ext Numeric vector with names "xmin", "xmax", "ymin", "ymax" in microns
#' indicating the spatial extent covered by the image. If \code{NULL}, then the
#' extent will be inferred from the metadata, from physical pixel size and the
#' number of pixels.
#' @param isFull Logical, if the extent specified in \code{ext} is the full
#' extent. If \code{ext = NULL} so it will be inferred from metadata then
#' \code{isFull = TRUE} will be set internally.
#' @param origin Origin of the image in the x-y plane, defaults to \code{c(0,0)}.
#' This is shifted when the image is translated.
#' @importFrom RBioFormats read.metadata coreMetadata
#' @return A \code{BioFormatsImage} object.
#' @name BioFormatsImage
#' @concept Image and raster
#' @seealso [isFull()], [origin()]
#' @export
setClass("BioFormatsImage", contains = "VirtualSpatialImage",
         slots = c(path = "character", ext = "numeric", isFull = "logical",
                   origin = "numeric"))

setValidity("BioFormatsImage", function(object) {
    outs <- character(3)
    e <- tryCatch(.check_bbox(object@ext), error = function(e) e$message)
    if (!is.null(e)) outs[1] <- e
    outs[2] <- if (is.na(object@isFull)) "isFull must be either TRUE or FALSE, not NA." else ""
    if (length(object@origin) != 2L || anyNA(object@origin) || !is.numeric(object@origin)) {
        outs[3] <- "origin must be a numeric vector of length 2 without NAs."
    }
    outs <- outs[outs != ""]
    if (length(outs)) return(outs) else TRUE
})

.get_fullres_scale_factor <- function(file) {
    check_installed(c("xml2", "RBioFormats"))
    xml_meta <- RBioFormats::read.omexml(file) |>
        xml2::read_xml() |> xml2::as_list()
    psx <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeX") |> as.numeric()
    psy <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeY") |> as.numeric()
    if (!length(psx) && !length(psy)) {
        stop("Physical pixel size absent from image metadata.")
    }
    #if (!is.null(psx) && is.null(psy)) {
    #    message("Assuming equal physical pixel size on x and y.")
    #    psy <- psx
    #} else if (is.null(psx) && !is.null(psy)) {
    #    message("Assuming equal physical pixel size on x and y.")
    #    psx <- psy
    #}
    sfx <- 1/psx
    sfy <- 1/psy
    c(sfx, sfy)
}

.get_fullres_size <- function(file) {
    #check_installed("RBioFormats")
    meta <- read.metadata(file) |>
        coreMetadata(series = 1L)
    sizeX_full <- meta$sizeX
    sizeY_full <- meta$sizeY
    c(sizeX_full, sizeY_full)
}

.get_full_ext <- function(file) {
    sfs <- .get_fullres_scale_factor(file)
    sfx <- sfs[1]; sfy <- sfs[2]
    size_full <- .get_fullres_size(file)
    sizeX_full <- size_full[1]; sizeY_full <- size_full[2]
    c(xmin = 0, ymin = 0, xmax = sizeX_full/sfx, ymax = sizeY_full/sfy)
}
# I think there should be separate documentation page for BioFormatsImage
#' @rdname BioFormatsImage
#' @export
BioFormatsImage <- function(path, ext = NULL, isFull = TRUE,
                            origin = c(0,0)) {
    # It would still be nice to automatically fill in ext from metadata
    # Previously I used NA to imply that it's full extent, but that could be confusing
    # I'll use another field, isFull, to indicate if it's full extent
    if (is.null(ext)) {
        ext <- .get_full_ext(path)
    }
    new("BioFormatsImage", path = path, ext = ext, isFull = isFull, origin = origin)
}

#' Other \code{BioFormatsImage} getters
#'
#' \code{isFULL} indicates if the extent is the full extent of the image.
#' \code{origin} gets the x-y coordinates of the origin of the image, i.e. the
#' smallest possible x-y coordinate values within the full image.
#'
#' @param x A \code{\link{BioFormatsImage}} object.
#' @return For \code{isFull}: Logical scalar indicating whether the extent is
#'   the full extent. For \code{origin}: Numeric vector of length 2.
#' @name BioFormatsImage-getters
#' @aliases isFull, origin
NULL

#' @rdname BioFormatsImage-getters
#' @export
setMethod("isFull", "BioFormatsImage", function(x) x@isFull)

#' @rdname BioFormatsImage-getters
#' @export
setMethod("origin", "BioFormatsImage", function(x) x@origin)

setReplaceMethod("origin", "BioFormatsImage", function(x, value) {
    x@origin <- value
    x
})

# EBImage==============

#' Representation of EBImage images in SFE objects
#'
#' This is a thin wrapper around the \code{\link{EBImage::Image}} class so it
#' inherits from \code{VirtualSpatialImage} to be compatible with
#' \code{SpatialExperiment} from which SFE inherits. An \code{ext} field is
#' added to specify the spatial extent of the image in microns to facilitate
#' geometric operations on the SFE object (including the images) and plotting
#' with \code{Voyager}.
#'
#' @inheritParams BioFormatsImage
#' @param img An \code{Image} object or anything that inherits from \code{Image}
#'   such as \code{\link{RBioFormats::AnnotatedImage}}.
#' @return An \code{EBImage} object.
#' @importClassesFrom EBImage Image
#' @name EBImage
#' @export
#' @concept Image and raster
setClass("EBImage", contains = "VirtualSpatialImage",
         slots = c(image = "Image", ext = "numeric"))

setValidity("EBImage", function(object) {
    out <- tryCatch(.check_bbox(object@ext), error = function(e) e$message)
    if (length(out)) return(out) else TRUE
})

#' @rdname EBImage
#' @export
EBImage <- function(img, ext = NULL) {
    if (is.null(ext)) stop("Extent must be specified for EBImage.")
    new("EBImage", image = img, ext = ext)
}

.shift_ext <- function(x, v) {
    x[c("xmin", "xmax")] <- x[c("xmin", "xmax")] + v[1]
    x[c("ymin", "ymax")] <- x[c("ymin", "ymax")] + v[2]
    x
}

# Coercion of Image classes===================
.toEBImage <- function(x, resolution = 4L) {
    # PhysicalSizeX seems to be a standard field
    if (length(resolution) != 1L ||
        !isTRUE(all.equal(floor(resolution), resolution))) {
        stop("resolution must be integer of length 1.")
    }
    file <- imgSource(x)
    sfs <- .get_fullres_scale_factor(file)
    sfx <- sfs[1]; sfy <- sfs[2]
    size_full <- .get_fullres_size(file)
    sizeX_full <- size_full[1]; sizeY_full <- size_full[2]
    min_nms <- c("xmin", "ymin")
    max_nms <- c("xmax", "ymax")
    x_nms <- c("xmin", "xmax")
    y_nms <- c("ymin", "ymax")
    if (!isFull(x)) {
        bbox_use <- ext(x) |> .shift_ext(v = -origin(x))
        # Convert to full res pixel space
        bbox_use[x_nms] <- bbox_use[x_nms] * sfx
        bbox_use[y_nms] <- bbox_use[y_nms] * sfy
        # Convert to lower res pixel space
        if (resolution > 1L) {
            meta <- read.metadata(file) |>
                coreMetadata(series = resolution)
            sfx2 <- meta$sizeX/sizeX_full
            sfy2 <- meta$sizeY/sizeY_full
            bbox_use[x_nms] <- bbox_use[x_nms] * sfx2
            bbox_use[y_nms] <- bbox_use[y_nms] * sfy2
        } else sfx2 <- sfy2 <- 1
        bbox_use[min_nms] <- floor(bbox_use[min_nms])
        bbox_use[max_nms] <- ceiling(bbox_use[max_nms])
        # For instance, say xmin = 0. Then it should start with pixel 1.
        subset_use <- list(x = seq(bbox_use["xmin"]+1L, bbox_use["xmax"]-1L, by = 1L),
                           y = seq(bbox_use["ymin"]+1L, bbox_use["ymax"]-1L, by = 1L))
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
    EBImage(img, .shift_ext(ext_use, origin(x)))
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
    x <- x@image
    if (dim(x)[3] == 3L) {
        names(x) <- c("r", "g", "b")
        # Remove RGB settings, better plot without it
        terra::RGB(x) <- NULL
    }
    x <- .resample_spat(x, maxcell)
    # Problem: in EBImage, the 1st dimension is x/column in the image
    if (dim(x)[3] == 3L)
        out <- terra::as.array(x) |> aperm(c(2,1,3)) |> Image(colormode = "Color")
    else
        out <- terra::as.array(x)[,,1] |> t() |> Image(colormode = "Grayscale")
    EBImage(out, ext(x) |> as.vector())
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
#' @seealso toSpatRasterImage
#' @aliases toEBImage
#' @export
#' @concept Image and raster
setMethod("toEBImage", "BioFormatsImage", .toEBImage)

#' @rdname toEBImage
#' @export
setMethod("toEBImage", "SpatRasterImage", .toEBImage2)

#' Convert images to SpatRasterImage
#'
#' The resolution specified from the OME-TIFF file will be read into memory and
#' written to disk as a GeoTIFF file that has the extent. The output file will
#' have the same file name as the input file except without the \code{ome} in
#' the extension.
#'
#' @inheritParams toEBImage
#' @param overwrite Logical, whether to overwrite existing file of the same
#'   name.
#' @param x Either a \code{BioFormatsImage} or \code{EBIImage} object.
#' @return A \code{SpatRasterImage} object
#' @aliases toSpatRasterImage
#' @name toSpatRasterImage
#' @seealso toEBImage
#' @export
#' @concept Image and raster
setMethod("toSpatRasterImage", "EBImage",
          function(x, file_out = "img.tiff", overwrite = FALSE) {
    m <- as.array(imgRaster(x))
    if (length(dim(m)) == 3L) m <- aperm(m, c(2,1,3))
    else m <- t(m)
    r <- rast(m, extent = ext(x))
    if (!file.exists(file_out) || overwrite) {
        message(">>> Saving image with `.tif` (non OME-TIFF) format:",
                paste0("\n", file_out))
        writeRaster(r, file_out, overwrite = overwrite)
    }
    rast(file_out) |> SpatRasterImage()
})

#' @rdname toSpatRasterImage
#' @export
setMethod("toSpatRasterImage", "BioFormatsImage",
          function(x, resolution = 4L, overwrite = FALSE) {
    #check_installed("RBioFormats")
    # Only for OME-TIFF, haven't tested on other BioFormats
    img <- toEBImage(x, resolution)
    img_fn <- gsub(".ome", paste0("_res", resolution), imgSource(x))
    toSpatRasterImage(img, img_fn, overwrite)
})

# Methods======================

#' Get and set extent of image objects
#'
#' Unlike in \code{SpatialExperiment}, images in SFE have extents which are used
#' to align them to the geometries and in geometric operations on SFE objects.
#' These functions get or set the extent for S4 image classes inheriting from
#' \code{VirtualSpatialImage} implemented in the SFE package.
#'
#' @param x A \code{*Image} object.
#' @param value A numeric vector with names "xmin", "xmax", "ymin", "ymax"
#' specifying the extent to use.
#' @return Getters return a numeric vector specifying the extent. Setters return
#' a \code{*Image} object of the same class as the input.
#' @name ext
#' @aliases ext
#' @concept Image and raster
NULL

#' @rdname ext
#' @export
setMethod("ext", "BioFormatsImage", function(x) x@ext[c("xmin", "xmax", "ymin", "ymax")])

#' @rdname ext
#' @export
setMethod("ext", "EBImage", function(x) x@ext[c("xmin", "xmax", "ymin", "ymax")])

#' @rdname ext
#' @export
setMethod("ext", "SpatRasterImage", function(x) ext(x@image) |> as.vector())

.set_ext <- function(x, value) {
    x@ext <- value[c("xmin", "xmax", "ymin", "ymax")]
    x
}

#' @rdname ext
#' @export
setReplaceMethod("ext", c("BioFormatsImage", "numeric"), .set_ext)


#' @rdname ext
#' @export
setReplaceMethod("ext", c("EBImage", "numeric"), .set_ext)

#' @rdname ext
#' @export
setReplaceMethod("ext", c("SpatRasterImage", "numeric"),
                 function(x, value) {
                     ext(x@image) <- value[c("xmin", "xmax", "ymin", "ymax")]
                     x
                 })


#' Methods for handling image-related data
#'
#' Generics of these functions are defined in \code{SpatialExperiment}, except
#' for \code{transposeImg}. These SFE methods cater to the new image-related
#' classes in SFE. The SPE method for \code{getImg}, \code{rmvImg}, and
#' \code{imgRaster} don't need to be modified for SFE and are hence not
#' implemented here, but are simply re-exported.
#'
#' Method of \code{\link{transposeImg}}, \code{\link{mirrorImg}}, and
#' \code{\link{rotateImg}} perform the method on all images within the SFE
#' object that are specified with \code{sample_id} and \code{image_id}.
#'
#' @inheritParams mirrorImg
#' @inheritParams rotateImg
#' @param x A SFE object.
#' @param sample_id Which sample the image is associated with. Use
#'   \code{\link{sampleIDs}} to get sample IDs present in the SFE object.
#' @param image_id Image ID, such as "lowres" and "hires" for Visium data and
#'   "DAPI" and "PolyT" for Vizgen MERFISH data.
#' @param file File from which to read the image.
#' @param extent A numeric vector of length 4 with names of the set xmin, ymin,
#'   xmax, and ymax, specifying the extent of the image.
#' @param scale_fct Scale factor -- multiply pixel coordinates in full
#'   resolution image by this scale factor should yield pixel coordinates in a
#'   different resolution. \code{extent} takes precedence over \code{scale_fct}.
#' @name SFE-image
#' @concept Image and raster
#' @family image methods
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' img_path <- system.file(file.path("extdata", "sample01", "outs", "spatial",
#' "tissue_lowres_image.png"), package = "SpatialFeatureExperiment")
#' sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres", scale_fct =
#' 0.023)
#' img <- getImg(sfe)
#' # SpatRasterImage method
#' img_t <- transposeImg(img)
#' # SFE method
#' sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "lowres")
NULL

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
        spi <- EBImage(img, extent)
    } else {
        if (!.path_valid2(img))
            stop("img is not a valid file path.")
        e <- tryCatch(suppressWarnings(rast(img)), error = function(e) e)
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
                # NOTE: potential bug when geometry is flipped
                img <- w
                flip <- FALSE
            }
            if (flip) img <- terra::flip(img)
            spi <- new("SpatRasterImage", image = img)
        }
    }

    DataFrame(
        sample_id,
        image_id,
        data=I(list(spi)),
        scaleFactor=scale_fct)
}

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
                  new <- lapply(old, transposeImg, resolution = resolution)
                  imgData(x)$data[idx] <- new
              }
              return(x)
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
                  new <- lapply(old, mirrorImg, resolution = resolution,
                                direction = direction)
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "SpatialFeatureExperiment",
          function(x, sample_id=NULL, image_id=NULL, degrees,
                   resolution = 4L, maxcell = 1e7) {
              sample_id <- .check_sample_id(x, sample_id, one = TRUE)
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  new <- lapply(old, rotateImg, resolution = resolution,
                                maxcell = maxcell,
                                degrees = degrees)
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

#' @rdname SFE-image
#' @export
setMethod("translateImg", "SpatialFeatureExperiment",
          function(x, sample_id=NULL, image_id=NULL, v) {
              sample_id <- .check_sample_id(x, sample_id, one = TRUE)
              old <- getImg(x, sample_id, image_id)
              if (!is.null(old)) {
                  if (!is.list(old)) old <- list(old)
                  idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
                  new <- lapply(old, translateImg, v = v)
                  imgData(x)$data[idx] <- new
              }
              return(x)
          })

# imgRaster----------

#' Get the image from *Image class
#'
#' In SFE, S4 classes inheriting from \code{VirtualSpatialImage} have been
#' implemented to make these image classes compatible with
#' \code{SpatialExperiment}. The \code{imgRaster} methods in SFE are meant to
#' extract the original image from the \code{*Image} classes, such as
#' \code{SpatRaster} from \code{SpatRasterImage}, and \code{Image} from
#' \code{EBImage} and \code{BioFormatsImage}. For \code{BioFormatsImage}, the
#' image of the specified resolution will be read into memory as
#' \code{AnnotatedImage}, which inherits from \code{EBImage::Image}.
#'
#' @param x An object of class \code{*Image} as implemented in this package.
#' @param resolution Resolution to read in from OME-TIFF, defaults to 4, which
#'   is a medium resolution in Xenium.
#' @return \code{SpatRaster} from \code{SpatRasterImage}, and \code{Image} from
#'   \code{EBImage} and \code{BioFormatsImage}. For \code{BioFormatsImage}, the
#'   image of the specified resolution will be read into memory as
#'   \code{AnnotatedImage}, which inherits from \code{EBImage::Image}.
#' @name imgRaster
#' @export
#' @concept Image and raster
#' @family image methods
setMethod("imgRaster", "SpatRasterImage", function(x) {
    if (is(x@image, "PackedSpatRaster")) unwrap(x@image)
    else x@image
})

#' @rdname imgRaster
#' @export
setMethod("imgRaster", "BioFormatsImage", function(x, resolution = 4L) {
    toEBImage(x, resolution) |> imgRaster()
})

#' @rdname imgRaster
#' @export
setMethod("imgRaster", "EBImage", function(x) x@image)

# TODO: imgRaster setter function since here I want to allow image processing
# like adjusting brightness and contrast, blurring, sharpening, opening, closing, and so on
# but what if it changes the extent?

# imgSource--------------

#' Source of images that are on disk
#'
#' Get the file path of images that are on disk and not read into memory. Only
#' applies to \code{SpatRasterImage} and \code{BioFormatsImage}.
#'
#' @inheritParams imgRaster
#' @return String, file path to the original image on disk. For
#'   \code{SpatRasterImage}, if the image is loaded into memory, then
#'   \code{NULL}.
#' @name imgSource
#' @export
#' @concept Image and raster
#' @importFrom terra sources
#' @family image methods
setMethod("imgSource",
          "SpatRasterImage",
          function(x) {
              sources(imgRaster(x))
          })

#' @rdname imgSource
#' @export
setMethod("imgSource", "BioFormatsImage", function(x) x@path)

setMethod("imgSource", "EBImage", function(x) NA_character_)

# Transpose-------------

#' Transpose images
#'
#' Swap rows and columns of images. In effect, this will flip the image around
#' the diagonal running from top left to bottom right.
#'
#' @inheritParams imgRaster
#' @param ... Ignored. It's there so different methods can all be passed to the
#' same \code{lapply} in the method for SFE objects.
#' @return For \code{SpatRasterImage} and \code{EBImage}, object of the same
#'   class. For \code{BioFormatsImage}, the image of the specified resolution is
#'   read into memory and then the \code{EBImage} method is called, returning
#'   \code{EBImage}. For the extent: xmin and xmax are switched with ymin and
#'   ymax.
#' @name transposeImg
#' @concept Image and raster
#' @export
#' @family image methods
setMethod("transposeImg", "SpatRasterImage",
          function(x, ...) {
              x@image <- terra::trans(imgRaster(x))
              # What terra does to extent: swap xmin and xmax with ymin and ymax
              x
          })

#' @rdname transposeImg
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

#' @rdname transposeImg
#' @export
setMethod("transposeImg", "EBImage",
          function(x, ...) {
              x@image <- EBImage::transpose(x@image)
              # Extent
              ext(x) <- .trans_extent(ext(x))
              x
          })

# Mirror-----------

#' Mirror/flip images
#'
#' Flip images along the middle horizontal or vertical axis.
#'
#' @inheritParams terra::flip
#' @inheritParams transposeImg
#' @return For \code{SpatRasterImage} and \code{EBImage}, object of the same
#'   class. For \code{BioFormatsImage}, the image of the specified resolution is
#'   read into memory and then the \code{EBImage} method is called, returning
#'   \code{EBImage}. The extent is not changed.
#' @name mirrorImg
#' @concept Image and raster
#' @export
#' @family image methods
setMethod("mirrorImg", "SpatRasterImage",
          function(x, direction = c("vertical", "horizontal"), ...) {
              direction <- match.arg(direction)
              x@image <- terra::flip(imgRaster(x), direction = direction)
              x
          })

#' @rdname mirrorImg
#' @export
setMethod("mirrorImg", "BioFormatsImage",
          function(x, direction = c("vertical", "horizontal"),
                   resolution = 4L) {
              direction <- match.arg(direction)
              x <- toEBImage(x, resolution)
              mirrorImg(x, direction)
          })

#' @rdname mirrorImg
#' @export
setMethod("mirrorImg", "EBImage",
          function(x, direction = c("vertical", "horizontal"), ...) {
              direction <- match.arg(direction)
              fun <- if (direction == "vertical") EBImage::flip else EBImage::flop
              x@image <- fun(x@image)
              x
          })

# Rotate---------------

# Just realized that SpatialExperiment only rotates by multiples of 90 degrees
# Can apply that and reconstruct rast. Will load into memory.

#' Rotate image
#'
#' As in \code{SpatialExperiment}, rotation here must be a multiple of 90
#' degrees.
#'
#' @inheritParams transposeImg
#' @param degrees How many degrees to rotate, must be multiples of 90. Positive
#'   number means clockwise and negative number means counterclockwise.
#' @param maxcell Max number of pixels to load \code{SpatRasterImage} into
#'   memory. The default 1e7 is chosen because this is the approximate number of
#'   pixels in the medium resolution image at \code{resolution = 4L} in Xenium
#'   OME-TIFF to make different methods of this function consistent.
#' @return An \code{EBImage} object. Both \code{SpatRasterImage} and
#'   \code{BioFormatsImage} will be loaded into memory as \code{EBImage} before
#'   rotating.
#' @name rotateImg
#' @concept Image and raster
#' @export
#' @family image methods
setMethod("rotateImg", "SpatRasterImage", # Deal with rotating SpatRaster?
          function(x, degrees, maxcell = 1e7, ...) {
              # Not sure what exactly to do. I think convert to EBImage as well.
              stopifnot(
                  length(degrees) == 1,
                  is.numeric(degrees),
                  degrees %% 90 == 0)
              x <- toEBImage(x, maxcell)
              rotateImg(x, degrees)
          })

#' @rdname rotateImg
#' @export
setMethod("rotateImg", "BioFormatsImage",
          function(x, degrees, resolution = 4L, ...) {
              # Only allow multiples of 90 degrees, deal with extent
              stopifnot(
                  length(degrees) == 1,
                  is.numeric(degrees),
                  degrees %% 90 == 0)
              x <- toEBImage(x, resolution)
              rotateImg(x, degrees)
          })

bbox_center <- function(bbox) {
    c(mean(bbox[c("xmin", "xmax")]),
      mean(bbox[c("ymin", "ymax")]))
}

#' @rdname rotateImg
#' @export
setMethod("rotateImg", "EBImage",
          function(x, degrees, ...) {
              stopifnot(
                  length(degrees) == 1,
                  is.numeric(degrees),
                  degrees %% 90 == 0)
              x@image <- EBImage::rotate(x@image, degrees)
              # Deal with extent
              if (degrees %% 180 > 0) {
                  ext_old <- ext(x)
                  center <- bbox_center(ext_old)
                  x_dist <- ext_old[["xmax"]] - center[1]
                  y_dist <- ext_old[["ymax"]] - center[2]
                  ext_new <- c(xmin = center[1] - y_dist,
                               xmax = center[1] + y_dist,
                               ymin = center[2] - x_dist,
                               ymax = center[2] + x_dist)
                  ext(x) <- ext_new
              }
              x
          })

# Translate-------------

#' Translate/shift image in space
#'
#' This function shifts the spatial extent of the image in the x-y plane.
#'
#' @inheritParams transposeImg
#' @param v Numeric vector of length 2 to shift the image in the x-y plane.
#' @return A \code{*Image} object of the same class that has been shifted in
#'   space.
#' @name translateImg
#' @aliases translateImg
#' @importFrom terra shift
#' @export
setMethod("translateImg", "SpatRasterImage", function(x, v) {
    img <- imgRaster(x)
    img <- shift(img, dx = v[1], dy = v[2])
    x@image <- img
    x
})

#' @rdname translateImg
#' @export
setMethod("translateImg", "BioFormatsImage", function(x, v) {
    ext(x) <- .shift_ext(ext(x), v)
    origin(x) <- origin(x) + v
    x
})

#' @rdname translateImg
#' @export
setMethod("translateImg", "EBImage", function(x, v) {
    ext(x) <- .shift_ext(ext(x), v)
    x
})

# Crop---------------

#' Crop images
#'
#' Crop images of class \code{*Image} in this package with a bounding box.
#'
#' @inheritParams transposeImg
#' @param bbox Numeric vector with names "xmin", "xmax", "ymin", "ymax", in any
#'   order, to specify the bounding box.
#' @return Image of the same class as input but cropped. For
#'   \code{BioFormatsImage}, the image is not loaded into memory; only the
#'   extent is changed.
#' @name cropImg
#' @aliases cropImg
#' @concept Image and raster
#' @export
#' @family image methods
setMethod("cropImg", "SpatRasterImage", function(x, bbox) {
    x@image <- terra::crop(x@image, bbox, snap = "out")
    x
})

#' @rdname cropImg
#' @export
setMethod("cropImg", "BioFormatsImage", function(x, bbox) {
    bbox_old <- ext(x)
    bbox_use <- c(xmin = max(bbox_old["xmin"], bbox["xmin"]),
                  xmax = min(bbox_old["xmax"], bbox["xmax"]),
                  ymin = max(bbox_old["ymin"], bbox["ymin"]),
                  ymax = min(bbox_old["ymax"], bbox["ymax"]))
    ext(x) <- bbox_use
    if (!isTRUE(all.equal(bbox_use, bbox_old[names(bbox_use)]))) {
        x@isFull <- FALSE
    }
    x
})

#' @rdname cropImg
#' @export
setMethod("cropImg", "EBImage", function(x, bbox) {
    # Convert bbox to pixel range based on ext(x)
    bbox_old <- ext(x)
    dim_old <- dim(x@image)
    sfx <- dim_old[1]/(bbox_old["xmax"] - bbox_old["xmin"])
    sfy <- dim_old[2]/(bbox_old["ymax"] - bbox_old["ymin"])
    bbox_use <- bbox
    bbox_use[c("xmin", "xmax")] <- bbox_use[c("xmin", "xmax")] * sfx
    bbox_use[c("ymin", "ymax")] <- bbox_use[c("ymin", "ymax")] * sfy
    bbox_use[c("xmin", "ymin")] <- floor(bbox_use[c("xmin", "ymin")])
    bbox_use[c("xmax", "ymax")] <- ceiling(bbox_use[c("xmax", "ymax")])

    if (length(dim_old) == 3L)
        x@image <- x@image[seq(bbox_use["xmin"]+1L, bbox_use["xmax"]-1L),
                           seq(bbox_use["ymin"]+1L, bbox_use["ymax"]-1L),]
    else
        x@image <- x@image[seq(bbox_use["xmin"]+1L, bbox_use["xmax"]-1L),
                           seq(bbox_use["ymin"]+1L, bbox_use["ymax"]-1L)]

    bbox_new <- bbox_use
    bbox_new[c("xmin", "xmax")] <- bbox_new[c("xmin", "xmax")] / sfx
    bbox_new[c("ymin", "ymax")] <- bbox_new[c("ymin", "ymax")] / sfy
    ext(x) <- bbox_new
    x
})

.crop_imgs <- function(x, bboxes) {
    # Crop all images across samples in an SFE object
    if (nrow(imgData(x))) {
        samples <- sort(sampleIDs(x))
        imgData(x) <- imgData(x)[order(imgData(x)$sample_id),]
        if (length(samples) == 1L) {
            bboxes <- matrix(bboxes, ncol = 1, dimnames = list(names(bboxes), samples))
        }
        new_imgs <- lapply(samples, function(s) {
            img_data <- imgData(x)$data[imgData(x)$sample_id == s]
            bbox_use <- bboxes[c("xmin", "xmax", "ymin", "ymax"),s]
            lapply(img_data, cropImg, bbox = bbox_use)
        })
        new_imgs <- unlist(new_imgs, recursive = FALSE)
        imgData(x)$data <- I(new_imgs)
    }
    x
}
