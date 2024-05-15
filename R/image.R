# Superclass for *Image with extent
#' @importFrom S4Vectors showAsCell
setClass("AlignedSpatialImage", contains = c("VirtualSpatialImage", "VIRTUAL"))
setMethod("NROW", "AlignedSpatialImage", function(x) 1L)
setMethod("NCOL", "AlignedSpatialImage", function(x) 1L)

# SpatRasterImage====================

#' SpatRaster representation of images in SFE objects
#'
#' \code{SpatialFeatureExperiment} and the \code{Voyager} package work with
#' images differently from \code{SpatialExperiment}. In SFE and
#' \code{Voyager}'s, plotting functions for SFE objects, the images can be read
#' with \code{\link{rast}} and represented as \code{SpatRaster}, so the image is
#' not entirely loaded into memory unless necessary. Plotting will not load a
#' large image into memory; rather the image will be downsampled and the
#' downsampled version is plotted. A \code{SpatRasterImage} object (as of Bioc
#' 3.19 or SFE version 1.6 and above) is a \code{SpatRaster} object but also
#' inheriting from \code{VirtualSpatialImage} as required by
#' \code{SpatialExperiment}.
#'
#' @param img A \code{\link{SpatRaster}} or \code{PackedSpatRaster} object.
#' @param object A \code{SpatRasterImage} object.
#' @return A \code{SpatRasterImage} object.
#' @importClassesFrom SpatialExperiment VirtualSpatialImage
#' @importFrom SpatialExperiment addImg mirrorImg imgData imgData<- imgRaster
#'   imgSource getImg rotateImg rmvImg
#' @importFrom terra ext ext<- inMemory
#' @importClassesFrom terra SpatRaster
#' @importClassesFrom EBImage Image
#' @importFrom methods setClassUnion
#' @concept Image classes
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

#' @rdname SpatRasterImage
#' @export
setClass("SpatRasterImage", contains=c("AlignedSpatialImage", "SpatRaster"))

#' @rdname SpatRasterImage
#' @export
SpatRasterImage <- function(img) new("SpatRasterImage", img)

setClass("PackedRasterImage", contains=c("AlignedSpatialImage", "PackedSpatRaster"))

#' @rdname SpatRasterImage
#' @export
setMethod("show", "SpatRasterImage", function(object) {
    d <- dim(object)
    dim <- paste(d[c(2,1,3)], collapse=" x ")
    str <- paste0(dim, " (width x height x channels) ", class(object), "\n")
    cat(str)
    str <- imgSource(object)
    if (!is.na(str)) {
        if (nchar(str) > 80) {
            str <- strwrap(str, width = 80)
        }
        cat("imgSource():\n ", str, "\n")
    }
})

setMethod("showAsCell", "SpatRasterImage", function(object) {
    d <- dim(object)
    paste(paste(d[c(2,1,3)], collapse=" x "), "SpatRasterImage")
})

# BioFormatsImage====================

#' On disk representation of BioFormats images in SFE object
#'
#' `r lifecycle::badge("experimental")` At present, the \code{BioFormatsImage}
#' is designed for OME-TIFF from Xenium and has not been tested on other formats
#' that can be read with \code{BioFormats}. The image is not loaded into memory,
#' and when it is, the the \code{BioFormatsImage} object is converted into
#' \code{\link{ExtImage}} because the loaded image is of a class that inherits
#' from \code{\link{Image}}. The \code{\link{ExtImage}} class is a thin wrapper
#' inheriting from \code{VirtualSpatialImage} so it's compatible with
#' \code{SpatialExperiment} from which SFE is derived. This class might
#' drastically change as it matures, say to accommodate other formats supported
#' by \code{BioFormats} and to store the transformation matrix rather than
#' loading image into memory upon transform.
#'
#' Spatial extent is inferred from OME-TIFF metadata if not specified. Physical
#' pixel size from the metadata is used to make the extent in micron space. If
#' physical pixel size is absent from metadata, then the extent will be in pixel
#' space, which might mean that the image will not align with the geometries
#' because often the geometry coordinates are in microns, so a warning is issued
#' in this case.
#'
#' Affine transformations can be specified in the \code{transformation}
#' argument, either by name or by directly specifying the matrix. The
#' transformations specified by name will always preserve the center of the
#' image. When named transformations are chained, name and parameter will be
#' converted to matrix and translation vector the second time a transformation
#' is specified. If the subsequent transformation happens to restore the image
#' to its original place, then transformation specifications will be removed.
#'
#' @param object A \code{BioFormatsImage} object.
#' @param path Path to an OME-TIFF image file.
#' @param ext Numeric vector with names "xmin", "xmax", "ymin", "ymax" in
#'   microns indicating the spatial extent covered by the image. If \code{NULL},
#'   then the extent will be inferred from the metadata, from physical pixel
#'   size and the number of pixels.
#' @param isFull Logical, if the extent specified in \code{ext} is the full
#'   extent. If \code{ext = NULL} so it will be inferred from metadata then
#'   \code{isFull = TRUE} will be set internally.
#' @param origin Origin of the whole image in the x-y plane, defaults to
#'   \code{c(0,0)}. This is shifted when the image is translated. This is not
#'   the same as xmin and xmax. For example, when the extent is only part of the
#'   whole image and the whole image itself can be spatially translated, the
#'   origin is needed to determine which part of the whole image this extent
#'   corresponds to.
#' @param transformation Named list specifying affine transformation. The list
#'   can have names "name" and named parameter of the transformation, e.g.
#'   \code{list(name = "mirror", direction = "vertical")}, "rotate" and degrees = 90
#'   (clockwise), and "scale" and factor = 2. The list can also have names "M"
#'   for a 2x2 linear transformation matrix in the xy plane and "v" for a
#'   translation vector of length 2 to specify general affine transformation.
#' @return A \code{BioFormatsImage} object.
#' @name BioFormatsImage
#' @aliases BioFormatsImage-class
#' @concept Image classes
#' @seealso [isFull()], [origin()]
#' @exportClass BioFormatsImage
setClass("BioFormatsImage", contains = "AlignedSpatialImage",
         slots = c(path = "character", ext = "numeric", isFull = "logical",
                   origin = "numeric", transformation = "list"))

.numeric2NA <- function(x) {
    length(x) != 2L || anyNA(x) || !is.numeric(x)
}
setValidity("BioFormatsImage", function(object) {
    outs <- character(4)
    e <- tryCatch(.check_bbox(object@ext), error = function(e) e$message)
    if (!is.null(e)) outs[1] <- e
    outs[2] <- if (is.na(object@isFull)) "isFull must be either TRUE or FALSE, not NA." else ""
    if (.numeric2NA(object@origin)) {
        outs[3] <- "origin must be a numeric vector of length 2 without NAs."
    }
    if (length(object@transformation)) {
        named <- "name" %in% names(object@transformation)
        mat <- setequal(names(object@transformation), c("M", "v"))
        if (!named && !mat) {
            outs[4] <- "transformation must be a list with names 'name' and 'param', or alternatively 'M' and 'v'."
        }
        if (named) {
            allowed_names <- c("transpose", "mirror", "rotate", "scale")
            if (!object@transformation$name %in% allowed_names) {
                outs[4] <- paste0("Name of transformation must be one of ",
                                  paste0(allowed_names, collapse = ", "))
                name <- object@transformation$name
            }
        } else if (mat) {
            M <- object@transformation$M
            v <- object@transformation$v
            if (!identical(dim(M), c(2,2)) || !is.numeric(M) || .numeric2NA(v)) {
                outs[4] <- paste0("M must be a 2x2 numeric matrix, and v a numeric vector of length 2 with no NA.")
            }
        }
    }
    outs <- outs[outs != ""]
    if (length(outs)) return(outs) else TRUE
})

.get_fullres_scale_factor <- function(file) {
    ps <- getPixelSize(file)
    psx <- ps[1]; psy <- ps[2]
    if (!length(ps)) {
        warning("Physical pixel size absent from image metadata, using pixel space.")
        sfx <- sfy <- 1
    } else {
        sfx <- 1/psx
        sfy <- 1/psy
    }
    c(sfx, sfy)
}

.get_fullres_size <- function(file) {
    check_installed("RBioFormats")
    coreMetadata <- RBioFormats::coreMetadata
    meta <- RBioFormats::read.metadata(file) |>
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

#' @rdname BioFormatsImage
#' @export
setMethod("show", "BioFormatsImage", function(object) {
    d <- dim(object)
    dim <- paste0("X: ", d[1], ", Y: ", d[2], ", C: ", d[3], ", Z: ", d[4], ", T: ", d[5])
    str <- paste0(dim, ", ", class(object), "\n")
    cat(str)
    str <- imgSource(object)
    if (!is.na(str)) {
        if (nchar(str) > 80) {
            str <- strwrap(str, width = 80)
        }
        cat("imgSource():\n ", str, "\n")
    }
})
setMethod("showAsCell", "BioFormatsImage", function(object) {
    paste(paste(dim(object)[1:3], collapse = " x "), "BioFormatsImage")
})

#' @rdname BioFormatsImage
#' @export
BioFormatsImage <- function(path, ext = NULL, isFull = TRUE, origin = c(0,0),
                            transformation = list()) {
    if (is.null(ext)) {
        if (isFull) ext <- .get_full_ext(path)
        else stop("Extent must be specified when isFull = FALSE")
    }
    new("BioFormatsImage", path = path, ext = ext, isFull = isFull,
        origin = origin, transformation = transformation)
}

# Combine transformations
.get_mv <- function(tr, bbox) {
    if ("name" %in% names(tr))
        do.call(.get_affine_mat, c(tr, list(bbox = bbox)))
    else tr
}
# The bbox argument here is used when the image is smaller than the overall
# bbox of the entire SFE object
.combine_transforms <- function(bfi, new, bbox = ext(bfi)) {
    old <- transformation(bfi)
    if (!length(old)) {
        nms <- c("xmin", "xmax", "ymin", "ymax")
        if (!isTRUE(all.equal(bbox[nms], ext(bfi)[nms]))) {
            new <- .get_mv(new, bbox)
        }
        transformation(bfi) <- new
    }else {
        old <- .get_mv(old, bbox)
        new <- .get_mv(new, bbox)
        out <- list(M = new$M %*% old$M, v = new$M %*% old$v + new$v)
        if (isTRUE(all.equal(out$M, diag(nrow = 2))) &&
            isTRUE(all.equal(as.vector(out$v), c(0,0)))) {
            out <- list()
        }
        transformation(bfi) <- out
    }
    bfi
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
#'   For \code{transformation}, a list.
#' @name BioFormatsImage-getters
#' @aliases isFull origin transformation
#' @concept Image classes
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

#' @rdname BioFormatsImage-getters
#' @export
setMethod("transformation", "BioFormatsImage", function(x) x@transformation)

setReplaceMethod("transformation", "BioFormatsImage", function(x, value) {
    x@transformation <- value
    x
})

# ExtImage==============

#' Use the EBImage \code{Image} class in SFE objects
#'
#' This is a thin wrapper around the \code{\link{Image}} class in the
#' \code{EBImage} package so it inherits from \code{VirtualSpatialImage} to be
#' compatible with \code{SpatialExperiment} from which SFE inherits. An
#' \code{ext} field is added to specify the spatial extent of the image in
#' microns to facilitate geometric operations on the SFE object (including the
#' images) and plotting with \code{Voyager}.
#'
#' @inheritParams BioFormatsImage
#' @param object An \code{ExtImage} object.
#' @param img An \code{Image} object or anything that inherits from \code{Image}
#'   such as \code{AnnotatedImage} in \code{RBioFormats}.
#' @return An \code{ExtImage} object.
#' @importClassesFrom EBImage Image
#' @importFrom EBImage Image
#' @name ExtImage
#' @aliases ExtImage-class
#' @exportClass ExtImage
#' @concept Image classes
setClass("ExtImage", contains = c("AlignedSpatialImage", "Image"),
         slots = c(ext = "numeric"))

setValidity("ExtImage", function(object) {
    out <- tryCatch(.check_bbox(object@ext), error = function(e) e$message)
    if (length(out)) return(out) else TRUE
})

#' @rdname ExtImage
#' @export
setMethod("show", "ExtImage", function(object) {
    d <- dim(object)
    dim <- paste(dim(object), collapse=" x ")
    if (length(d) == 2L)
        str <- paste0(dim, " (width x height) ", class(object), "\n")
    else if (length(d) == 3L)
        str <- paste0(dim, " (width x height x channels) ", class(object), "\n")
    cat(str)
})

setMethod("showAsCell", "ExtImage", function(object) {
    dim <- paste(dim(object), collapse=" x ")
    paste(dim, "ExtImage")
})

#' @rdname ExtImage
#' @export
ExtImage <- function(img, ext = NULL) {
    if (is.null(ext)) stop("Extent must be specified for ExtImage.")
    new("ExtImage", img, ext = ext)
}

.shift_ext <- function(x, v) {
    x[c("xmin", "xmax")] <- x[c("xmin", "xmax")] + v[1]
    x[c("ymin", "ymax")] <- x[c("ymin", "ymax")] + v[2]
    x
}

# Coercion of Image classes===================
.get_subset <- function(bbox_use, sizeY, channel = NULL) {
    ymin_px <- sizeY - bbox_use["ymax"]
    ymax_px <- sizeY - bbox_use["ymin"]
    bbox_img <- bbox_use
    bbox_img["ymin"] <- ymin_px; bbox_img["ymax"] <- ymax_px
    min_nms <- c("xmin", "ymin")
    max_nms <- c("xmax", "ymax")
    bbox_img[min_nms] <- floor(bbox_img[min_nms])
    bbox_img[max_nms] <- ceiling(bbox_img[max_nms])
    # For instance, say xmin = 0. Then it should start with pixel 1.
    list(x = seq(bbox_img["xmin"]+1L, bbox_img["xmax"]-1L, by = 1L),
         y = seq(bbox_img["ymin"]+1L, bbox_img["ymax"]-1L, by = 1L),
         c = channel)
}

.toExtImage <- function(x, resolution = 4L, channel = NULL) {
    check_installed("RBioFormats")
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

    coreMetadata <- RBioFormats::coreMetadata
    m <- RBioFormats::read.metadata(file)
    cm <- coreMetadata(m)
    if ("sizeX" %in% names(cm)) {
        # Indicating there's only 1 series/resolution
        resolution <- 1L
    } else {
        n_res <- RBioFormats::seriesCount(m)
        if (resolution > n_res) {
            message("Resolution subscript out of bound, reading the lowest resolution")
            resolution <- n_res
        }
        meta <- coreMetadata(m, series = resolution)
    }
    # Extent of lower resolution may not be the same as the top one
    # Need to more accurately infer extent
    # Infer the scale factor, and then get the difference from the rounding
    if (resolution > 1L) {
        fct_x <- sizeX_full/meta$sizeX
        fct_y <- sizeY_full/meta$sizeY
        fct_round <- round(fct_x) # Should be the same for x and y
        fctx2 <- fct_x/fct_round
        fcty2 <- fct_y/fct_round
        # The shift is worse when approaching the lower right corner
        o <- origin(x)
        o[1] <- o[1]/fctx2
        o[2] <- o[2]/fcty2
        origin(x) <- o
    } else fctx2 <- fcty2 <- 1

    if (!isFull(x)) {
        bbox_use <- .ext_(x) |> .shift_ext(v = -origin(x))
        # Convert to full res pixel space
        bbox_use[x_nms] <- bbox_use[x_nms] * sfx
        bbox_use[y_nms] <- bbox_use[y_nms] * sfy
        # Convert to lower res pixel space
        if (resolution > 1L) {
            sfx2 <- meta$sizeX*fctx2/sizeX_full
            sfy2 <- meta$sizeY*fcty2/sizeY_full
            bbox_use[x_nms] <- bbox_use[x_nms] * sfx2
            bbox_use[y_nms] <- bbox_use[y_nms] * sfy2
        } else sfx2 <- sfy2 <- 1
        # ext(x) has origin at the bottom left, while image indices have
        # origin at the top left. Need to convert the y indices.
        subset_use <- .get_subset(bbox_use, meta$sizeY, channel)
        bbox_use[min_nms] <- floor(bbox_use[min_nms])
        bbox_use[max_nms] <- ceiling(bbox_use[max_nms])
        # Extent should account for the pixels
        ext_use <- bbox_use
        ext_use[x_nms] <- ext_use[x_nms]/(sfx*sfx2)
        ext_use[y_nms] <- ext_use[y_nms]/(sfy*sfy2)
    } else {
        ext_use <- c(xmin = 0, ymin = 0, xmax = sizeX_full/(sfx*fctx2),
                     ymax = sizeY_full/(sfy*fcty2))
        subset_use <- list(c = channel)
    }
    img <- RBioFormats::read.image(file = file,
                                   resolution = resolution,
                                   filter.metadata = TRUE,
                                   read.metadata = FALSE,
                                   normalize = FALSE,
                                   subset = subset_use)
    out <- ExtImage(img, .shift_ext(ext_use, origin(x)))
    if (length(transformation(x))) {
        tr <- transformation(x)
        if (!"name" %in% names(tr)) tr$name <- "affine"
        trans_fun <- .get_img_fun(tr$name)
        name_param <- setdiff(names(tr), "name")
        out <- do.call(trans_fun, c(list(x = out), list(bbox_all = ext(out)), tr[name_param]))
    }
    out
}

# SpatRasterImage method: allow downsampling first
# From tidyterra and ggspavis::plotVisium
.resample_spat <- function(r, maxcell = 50000) {
    if (is.null(maxcell)) return(r)
    if (terra::ncell(r) > 1.1 * maxcell) {
        r <- terra::spatSample(r, maxcell,
                               as.raster = TRUE,
                               method = "regular"
        )
    }
    return(r)
}

.toExtImage2 <- function(x, maxcell = 1e7, channel = NULL) {
    # 1e7 comes from the number of pixels in resolution = 4L in the ome.tiff
    if (dim(x)[3] == 3L) {
        names(x) <- c("r", "g", "b")
        # Remove RGB settings, better plot without it
        terra::RGB(x) <- NULL
    }
    if (is.null(maxcell) && !inMemory(x)) {
        stop("Image cannot be loaded into memory. Please downsample by setting maxcell.")
    }
    x <- .resample_spat(x, maxcell)
    # Problem: in ExtImage, the 1st dimension is x/column in the image
    if (dim(x)[3] == 3L)
        out <- terra::as.array(x) |> aperm(c(2,1,3)) |> Image(colormode = "Color")
    else
        out <- terra::as.array(x)[,,1] |> t() |> Image(colormode = "Grayscale")
    ExtImage(out, ext(x))
}

#' Convert images to ExtImage
#'
#' The \code{ExtImage} class is a thin wrapper around the \code{Image} class in
#' \code{ExtImage} so it inherits from \code{VirtualSpatialImage} as required by
#' \code{SpatialExperiment} and has extent as used in Voyager's plotting
#' functions. This function converts \code{SpatRasterImage} (thin wrapper around
#' the class in \code{terra}) and \code{BioFormatsImage} into \code{ExtImage} for
#' image operations as implemented in the \code{ExtImage} package.
#'
#' @param x Either a \code{BioFormatsImage} or \code{SpatRasterImage} object.
#' @param resolution Integer, which resolution in the \code{BioFormatsImage} to
#'   read and convert. Defaults to 4, which is a lower resolution. Ignored if
#'   only 1 resolution is present.
#' @param channel Integer vector to indicate channel(s) to read. If \code{NULL},
#'   then all channels will be read.
#' @param maxcell Maximum number of pixels when \code{SpatRasterImage} is read
#'   into memory.
#' @return A \code{ExtImage} object. The image is loaded into memory.
#' @name toExtImage
#' @seealso toSpatRasterImage
#' @aliases toExtImage
#' @export
#' @concept Image classes
NULL

#' @rdname toExtImage
#' @export
setMethod("toExtImage", "BioFormatsImage", .toExtImage)

#' @rdname toExtImage
#' @export
setMethod("toExtImage", "SpatRasterImage", .toExtImage2)

#' Convert images to SpatRasterImage
#'
#' The resolution specified from the OME-TIFF file will be read into memory and
#' written to disk as a GeoTIFF file that has the extent. The output file will
#' have the same file name as the input file except without the \code{ome} in
#' the extension.
#'
#' @inheritParams toExtImage
#' @param overwrite Logical, whether to overwrite existing file of the same
#'   name.
#' @param save_geotiff Logical, whether to save the image to GeoTIFF file.
#' @param x Either a \code{BioFormatsImage} or \code{EBIImage} object.
#' @param file_out File to save the non-OME TIFF file for \code{SpatRaster}.
#' @return A \code{SpatRasterImage} object
#' @aliases toSpatRasterImage
#' @name toSpatRasterImage
#' @seealso toExtImage
#' @export
#' @concept Image classes
NULL

#' @rdname toSpatRasterImage
#' @export
setMethod("toSpatRasterImage", "ExtImage",
          function(x, save_geotiff = TRUE, file_out = "img.tiff", overwrite = FALSE) {
    m <- as.array(imgRaster(x))
    if (length(dim(m)) == 3L) m <- aperm(m, c(2,1,3))
    else m <- t(m)
    r <- rast(m, extent = ext(x))
    if (!save_geotiff) return(SpatRasterImage(r))
    if (!file.exists(file_out) || overwrite) {
        message(">>> Saving image with `.tiff` (non OME-TIFF) format:",
                paste0("\n", file_out))
        terra::writeRaster(r, file_out, overwrite = overwrite)
    }
    rast(file_out) |> SpatRasterImage()
})

#' @rdname toSpatRasterImage
#' @export
setMethod("toSpatRasterImage", "BioFormatsImage",
          function(x, save_geotiff = TRUE, resolution = 4L, channel = NULL,
                   overwrite = FALSE) {
    #check_installed("RBioFormats")
    # Only for OME-TIFF, haven't tested on other BioFormats
    img <- toExtImage(x, resolution)
    img_fn <- gsub("\\.(ome\\.)?tiff?$", paste0("_res", resolution,".tiff"), imgSource(x))
    toSpatRasterImage(img, save_geotiff, img_fn, overwrite)
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
#'   specifying the extent to use.
#' @note For \code{SpatRasterImage}, the image may be may not be loaded into
#' memory. You can check if the image is loaded into memory with
#' \code{terra::inMemory(imgRaster(x))}, and check the original file path with
#' \code{\link{imgSource}}. If the image is not loaded into memory, then the
#' original file must be present at the path indicated by
#' \code{\link{imgSource}} in order for any code using the image to work, which
#' includes this function \code{ext}.
#' @return Getters return a numeric vector specifying the extent. Setters return
#'   a \code{*Image} object of the same class as the input.
#' @note
#' For \code{BioFormatsImage}, internally only the pre-transform extent is
#' stored. The \code{ext} getter will apply the transformation on the fly. The
#' setter sets the pre-transformation extent.
#' @name ext
#' @aliases ext
#' @concept Image methods
#' @family image methods
NULL

.ext_ <- function(x) x@ext[c("xmin", "xmax", "ymin", "ymax")]

#' @rdname ext
#' @export
setMethod("ext", "BioFormatsImage",
          function(x) {
              out <- .transform_bbox(.ext_(x), transformation(x))
              out[c("xmin", "xmax", "ymin", "ymax")]
          })

#' @rdname ext
#' @export
setMethod("ext", "ExtImage", function(x) x@ext[c("xmin", "xmax", "ymin", "ymax")])

#' @rdname ext
#' @export
setMethod("ext", "SpatRasterImage", function(x) ext(as(x, "SpatRaster")) |> as.vector())

.set_ext <- function(x, value) {
    x@ext <- value[c("xmin", "xmax", "ymin", "ymax")]
    x
}

#' @rdname ext
#' @export
setReplaceMethod("ext", c("BioFormatsImage", "numeric"), .set_ext)


#' @rdname ext
#' @export
setReplaceMethod("ext", c("ExtImage", "numeric"), .set_ext)

#' @rdname ext
#' @export
setReplaceMethod("ext", c("SpatRasterImage", "numeric"),
                 function(x, value) {
                     x <- as(x, "SpatRaster")
                     ext(x) <- value[c("xmin", "xmax", "ymin", "ymax")]
                     SpatRasterImage(x)
                 })
# dim - BioFormatsImage----------

#' Find dimension of BioFormatsImage
#'
#' This is different from other classes. The metadata is read where the
#' dimensions in pixels can be found. The image itself is not read into memory
#' here.
#'
#' @param x A \code{\link{BioFormatsImage}} object.
#' @return An integer vector of length 5 showing the number of rows and columns
#'   in the full resolution image. The 5 dimensions are in the order of XYCZT:
#'   x, y, channel, z, and time. This is not changed by transformations. Use
#'   \code{\link{ext}} to see the extent after transformation.
#' @concept Image methods
#' @family image methods
#' @export
setMethod("dim", "BioFormatsImage", function(x) {
    check_installed("RBioFormats")
    coreMetadata <- RBioFormats::coreMetadata
    meta <- RBioFormats::read.metadata(imgSource(x)) |>
        coreMetadata(series = 1L)
    c(X=meta$sizeX, Y=meta$sizeY, C=meta$sizeC, Z=meta$sizeZ, "T"=meta$sizeT)
})

#' Image setter
#'
#' Modify or replace images stored in a \code{SpatialExperiment} object. This is
#' different from \code{\link{addImg}} which adds the image from files and can't
#' replace existing images, which is there to be consistent with
#' \code{SpatialExperiment}. This setter here can replace existing images with
#' another object that inherits from \code{VirtualSpatialImage}, including
#' \code{\link{SpatRasterImage}}, \code{\link{BioFormatsImage}}, and
#' \code{\link{ExtImage}}.
#'
#' @inheritParams SFE-image
#' @param x A \code{SpatialExperiment} object, which includes SFE.
#' @param scale_fct Scale factor to convert pixels in lower resolution to those
#'   in the full resolution. Only relevant to image classes implemented in
#'   \code{SpatialExperiment} but not \code{SpatialFeatureExperiment} because
#'   the spatial extent of images in SFE takes precedence.
#' @param value New version of image to add, must inherit from
#'   \code{VirtualSpatialImage}.
#' @return SFE object with the new image added.
#' @concept Image methods
#' @export
#' @importFrom methods signature
#' @importFrom SpatialExperiment imgData<-
#' @aliases Img<-
#' @examples
#' library(EBImage)
#' library(SFEData)
#' library(RBioFormats)
#' fp <- tempdir()
#' fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
#' # Weirdly the first time I get the null pointer error
#' try(sfe <- readXenium(fn))
#' sfe <- readXenium(fn)
#' img <- getImg(sfe) |> toExtImage(resolution = 1L)
#' img <- img[,,1] > 500
#' Img(sfe, image_id = "mask") <- img
#' imageIDs(sfe)
#' unlink(fn, recursive = TRUE)
#'
setMethod("Img<-", signature = "SpatialExperiment",
          function(x, sample_id = 1L, image_id, scale_fct = 1, value) {
              sample_id <- .check_sample_id(x, sample_id)
              df <- imgData(x)
              ind <- which(df$sample_id == sample_id & df$image_id == image_id)
              if (length(ind)) {
                  df$data[[ind]] <- value
              } else {
                  df_new <- DataFrame(
                      sample_id,
                      image_id,
                      data=I(list(value)),
                      scaleFactor=scale_fct)
                  df <- rbind(df, df_new)
              }
              imgData(x) <- df
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
#' object that are specified with \code{sample_id} and \code{image_id}. For
#' images that are not loaded into memory, \code{rotateImg} will load
#' \code{SpatRasterImage} into memory and all image operations except translate
#' will load \code{BioFormatsImage} into memory.
#'
#' @inheritParams mirrorImg
#' @inheritParams rotateImg
#' @inheritParams scaleImg
#' @inheritParams affineImg
#' @inheritParams SpatialExperiment::addImg
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
#' @param v A numeric vector of length 2 specifying the vector in the xy plane
#'   to translate the SFE object.
#' @note If the image is already a GeoTIFF file that already has an extent, then
#' the extent associated with the file will be honored and the \code{extent} and
#' \code{scale_fct} arguments are ignored. Transposing the image is just like
#' transposing a matrix. It's flipped about the line going from the top left to
#' the bottom right.
#' @name SFE-image
#' @concept Image methods
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
          function(x, imageSource, sample_id = 1L, image_id,
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
    # ExtImage
    if (is(img, "Image")) {
        if (is(img, "ExtImage")) spi <- img else spi <- ExtImage(img, extent)
    } else {
        if (!.path_valid2(img))
            stop("img is not a valid file path.")
        e <- tryCatch(suppressWarnings(rast(img)), error = function(e) e)
        if (is(e, "error") || grepl("\\.ome\\.tif", img)) {
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
            spi <- new("SpatRasterImage", img)
            if (flip) spi <- mirrorImg(spi, direction = "vertical")
        }
    }

    DataFrame(
        sample_id,
        image_id,
        data=I(list(spi)),
        scaleFactor=scale_fct)
}

.transform_img_sfe <- function(img_fun) {
    function(x, sample_id = 1L, image_id = NULL, ...) {
        sample_id <- .check_sample_id(x, sample_id, one = TRUE)
        old <- getImg(x, sample_id, image_id)
        if (!is.null(old)) {
            if (!is.list(old)) old <- list(old)
            idx <- SpatialExperiment:::.get_img_idx(x, sample_id, image_id)
            new <- lapply(old, img_fun, ...)
            imgData(x)$data[idx] <- new
        }
        return(x)
    }
}

#' @rdname SFE-image
#' @export
setMethod("transposeImg", "SpatialFeatureExperiment",
          function(x, sample_id = 1L, image_id = NULL,
                   maxcell = 1e7, filename = "") {
              .transform_img_sfe(transposeImg)(x, sample_id, image_id,
                                               maxcell = maxcell,
                                               filename = filename)
          })

#' @rdname SFE-image
#' @export
setMethod("mirrorImg", "SpatialFeatureExperiment",
          function(x, sample_id = 1L, image_id = NULL, direction = "vertical",
                   maxcell = 1e7, filename = "") {
              .transform_img_sfe(mirrorImg)(x, sample_id, image_id,
                                            maxcell = maxcell,
                                            filename = filename,
                                            direction = direction)
          })

#' @rdname SFE-image
#' @export
setMethod("rotateImg", "SpatialFeatureExperiment",
          function(x, sample_id = 1L, image_id = NULL, degrees,
                   maxcell = 1e7) {
              .transform_img_sfe(rotateImg)(x, sample_id, image_id,
                                            maxcell = maxcell,
                                            degrees = degrees)
          })

#' @rdname SFE-image
#' @export
setMethod("translateImg", "SpatialFeatureExperiment",
          function(x, sample_id = 1L, image_id = NULL, v) {
              .transform_img_sfe(translateImg)(x, sample_id, image_id,
                                               v = v)
          })

#' @rdname SFE-image
#' @export
setMethod("scaleImg", "SpatialFeatureExperiment",
          function(x, sample_id = 1L, image_id = NULL, factor) {
              .transform_img_sfe(scaleImg)(x, sample_id, image_id,
                                           factor = factor)
          })

#' @rdname SFE-image
#' @export
setMethod("affineImg", "SpatialFeatureExperiment",
          function(x, sample_id = 1L, image_id = NULL, M, v) {
              .transform_img_sfe(scaleImg)(x, sample_id, image_id,
                                           M = M, v = v)
          })

# imgRaster----------

#' Get the image from *Image class
#'
#' In SFE, S4 classes inheriting from \code{VirtualSpatialImage} have been
#' implemented to make these image classes compatible with
#' \code{SpatialExperiment}. The \code{imgRaster} methods in SFE are meant to
#' extract the original image from the \code{*Image} classes, such as
#' \code{SpatRaster} from \code{SpatRasterImage}, and \code{Image} from
#' \code{ExtImage} and \code{BioFormatsImage}. For \code{BioFormatsImage}, the
#' image of the specified resolution will be read into memory as
#' \code{AnnotatedImage}, which inherits from \code{EBImage::Image}.
#'
#' @param x An object of class \code{*Image} as implemented in this package.
#' @param resolution Resolution to read in from OME-TIFF, defaults to 4, which
#'   is a medium resolution in Xenium.
#' @return \code{SpatRaster} from \code{SpatRasterImage}, and \code{Image} from
#'   \code{ExtImage} and \code{BioFormatsImage}. For \code{BioFormatsImage}, the
#'   image of the specified resolution will be read into memory as
#'   \code{AnnotatedImage} and \code{ExtImage}, which both inherit from
#'   \code{EBImage::Image}.
#' @export
#' @name imgRaster
#' @aliases imgRaster,SpatRasterImage-method imgRaster,BioFormatsImage-method
#'   imgRaster,ExtImage-method
#' @concept Image methods
#' @family image methods
NULL

#' @export
setMethod("imgRaster", "SpatRasterImage", function(x) as(x, "SpatRaster"))

#' @export
setMethod("imgRaster", "BioFormatsImage", function(x, resolution = 4L) {
    toExtImage(x, resolution) |> imgRaster()
})

#' @export
setMethod("imgRaster", "ExtImage", function(x) as(x, "Image"))

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
#' @concept Image methods
#' @importFrom terra sources
#' @family image methods
NULL

#' @rdname imgSource
#' @export
setMethod("imgSource",
          "SpatRasterImage",
          function(x) {
              out <- sources(imgRaster(x))
              if (out == "") out <- NA_character_
              return(out)
          })

#' @rdname imgSource
#' @export
setMethod("imgSource", "BioFormatsImage", function(x) x@path)

#' @rdname imgSource
#' @export
setMethod("imgSource", "ExtImage", function(x) NA_character_)

# Transpose-------------

#' Transpose images
#'
#' Swap rows and columns of images. In effect, this will flip the image around
#' the diagonal running from top left to bottom right.
#'
#' @inheritParams imgRaster
#' @param filename Output file name for transformed SpatRaster.
#' @param maxcell Max number of pixels to load \code{SpatRasterImage} into
#'   memory. The default 1e7 is chosen because this is the approximate number of
#'   pixels in the medium resolution image at \code{resolution = 4L} in Xenium
#'   OME-TIFF to make different methods of this function consistent.
#' @param ... Ignored. It's there so different methods can all be passed to the
#'   same \code{lapply} in the method for SFE objects. Some methods have extra
#'   arguments.
#' @return For \code{SpatRasterImage} and \code{ExtImage}, object of the same
#'   class. For \code{BioFormatsImage}, the image of the specified resolution is
#'   read into memory and then the \code{ExtImage} method is called, returning
#'   \code{ExtImage}. For the extent: xmin and xmax are switched with ymin and
#'   ymax.
#' @name transposeImg
#' @aliases transposeImg
#' @concept Image affine transformation
#' @export
#' @family image methods
NULL

# From tidyterra
.resample_spat <- function(r, maxcell) {
    if (terra::ncell(r) > 1.1 * maxcell) {
        r <- terra::spatSample(r, maxcell,
                               as.raster = TRUE,
                               method = "regular"
        )
    }
    return(r)
}

setMethod(".transpose_img", "SpatRasterImage",
          function(x, bbox_all = ext(x), filename = "", maxcell = NULL) {
              if (!is.null(maxcell)) x <- .resample_spat(x, maxcell)
              ext_orig <- ext(x)
              x <- terra::trans(x, filename = filename)
              # Shift extent for overall bbox
              ext(x) <- .transform_bbox(ext_orig, tr = list(name = "transpose"),
                                        bbox_all = bbox_all)
              x
          })

#' @rdname transposeImg
#' @export
setMethod("transposeImg", "SpatRasterImage",
          function(x, filename = "", maxcell = NULL, ...) {
              .transpose_img(x, filename = filename, maxcell = maxcell,
                             bbox_all = ext(x))
          })

setMethod(".transpose_img", "BioFormatsImage",
          function(x, bbox_all = ext(x), ...) {
              .combine_transforms(x, list(name = "transpose"), bbox = bbox_all)
          })

#' @rdname transposeImg
#' @export
setMethod("transposeImg", "BioFormatsImage",
          function(x, ...) {
              .transpose_img(x, bbox_all = ext(x))
          })

setMethod(".transpose_img", "ExtImage",
          function(x, bbox_all = ext(x), ...) {
              x <- EBImage::transpose(x)
              # Extent
              ext(x) <- .transform_bbox(ext(x), list(name="transpose"),
                                        bbox_all = bbox_all)
              x
          })

#' @rdname transposeImg
#' @export
setMethod("transposeImg", "ExtImage",
          function(x, ...) {
              .transpose_img(x, bbox_all = ext(x))
          })

# Mirror-----------

#' Mirror/flip images
#'
#' Flip images along the middle horizontal or vertical axis.
#'
#' @inheritParams terra::flip
#' @inheritParams transposeImg
#' @return \code{*Image} object of the same class.
#' @name mirrorImg
#' @aliases mirrorImg
#' @concept Image affine transformation
#' @export
#' @family image methods
NULL

setMethod(".mirror_img", "SpatRasterImage",
          function(x, direction = c("vertical", "horizontal"),
                   bbox_all = NULL, filename = "", maxcell = NULL) {
              direction <- match.arg(direction)
              if (!is.null(maxcell)) x <- .resample_spat(x, maxcell)
              x <- terra::flip(imgRaster(x), direction = direction,
                                     filename = filename) |> SpatRasterImage()
              # Shift extent for overall bbox
              if (!is.null(bbox_all)) {
                  ext(x) <- .transform_bbox(ext(x),
                                            tr = list(name = "mirror",
                                                      direction = direction),
                                            bbox_all = bbox_all)
              }
              x
          })

#' @rdname mirrorImg
#' @export
setMethod("mirrorImg", "SpatRasterImage",
          function(x, direction = c("vertical", "horizontal"), filename = "",
                   maxcell = NULL, ...) {
              .mirror_img(x, direction = direction, filename = filename,
                          maxcell = maxcell, bbox_all = NULL)
          })

setMethod(".mirror_img", "BioFormatsImage",
          function(x, direction = c("vertical", "horizontal"),
                   bbox_all = ext(x), ...) {
              direction <- match.arg(direction)
              .combine_transforms(x, list(name = "mirror", direction = direction),
                                  bbox = bbox_all)
          })

#' @rdname mirrorImg
#' @export
setMethod("mirrorImg", "BioFormatsImage",
          function(x, direction = c("vertical", "horizontal"), ...) {
              .mirror_img(x, direction = direction, bbox_all = ext(x))
          })

setMethod(".mirror_img", "ExtImage",
          function(x, direction = c("vertical", "horizontal"),
                   bbox_all = NULL, ...) {
              direction <- match.arg(direction)
              fun <- if (direction == "vertical") EBImage::flip else EBImage::flop
              x <- fun(x)
              # Shift extent for overall bbox
              if (!is.null(bbox_all)) {
                  ext(x) <- .transform_bbox(ext(x),
                                            tr = list(name = "mirror",
                                                      direction = direction),
                                            bbox_all = bbox_all)
              }
              x
          })

#' @rdname mirrorImg
#' @export
setMethod("mirrorImg", "ExtImage",
          function(x, direction = c("vertical", "horizontal"), ...) {
              .mirror_img(x, direction = direction, bbox_all = NULL)
          })

# Rotate---------------

# Just realized that SpatialExperiment only rotates by multiples of 90 degrees
# Can apply that and reconstruct rast. Will load \code{SpatRasterImage} into
# memory.

#' Rotate image
#'
#' As in \code{SpatialExperiment}, rotation here must be a multiple of 90
#' degrees.
#'
#' @inheritParams transposeImg
#' @param degrees How many degrees to rotate. Positive number means clockwise
#'   and negative number means counterclockwise.
#' @return \code{SpatRasterImage} will be loaded into memory and converted to
#'   \code{ExtImage}. Otherwise \code{*Image} object of the same class.
#' @name rotateImg
#' @aliases rotateImg
#' @concept Image affine transformation
#' @export
#' @family image methods
NULL

setMethod(".rotate_img", "SpatRasterImage",
          function(x, degrees, bbox_all = ext(x), maxcell = 1e7) {
              # Not sure what exactly to do. I think convert to ExtImage as well.
              # Covert to BioFormatsImage if it's not in memory
              x <- toExtImage(x, maxcell)
              .rotate_img(x, degrees, bbox_all)
          })

#' @rdname rotateImg
#' @export
setMethod("rotateImg", "SpatRasterImage",
          function(x, degrees, maxcell = 1e7, ...) {
              .rotate_img(x, degrees, maxcell = maxcell, bbox_all = ext(x))
          })

setMethod(".rotate_img", "BioFormatsImage",
          function(x, degrees, bbox_all = ext(x), ...) {
              .combine_transforms(x, list(name = "rotate", degrees = degrees),
                                  bbox = bbox_all)
          })

#' @rdname rotateImg
#' @export
setMethod("rotateImg", "BioFormatsImage",
          function(x, degrees, ...) {
              .rotate_img(x, degrees = degrees, bbox_all = ext(x))
          })

setMethod(".rotate_img", "ExtImage",
          function(x, degrees, bbox_all = ext(x), ...) {
              x <- EBImage::rotate(x, degrees)
              ext(x) <- .transform_bbox(ext(x), list(name="rotate", degrees=degrees),
                                        bbox_all = bbox_all)
              x
          })

#' @rdname rotateImg
#' @export
setMethod("rotateImg", "ExtImage",
          function(x, degrees, ...) {
              .rotate_img(x, degrees = degrees, bbox_all = ext(x))
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
#' @concept Image affine transformation
#' @family image methods
#' @export
NULL

#' @rdname translateImg
#' @export
setMethod("translateImg", "SpatRasterImage", function(x, v, ...) {
    img <- imgRaster(x)
    img <- shift(img, dx = v[1], dy = v[2]) |> SpatRasterImage()
    x <- img
    x
})

#' @rdname translateImg
#' @export
setMethod("translateImg", "BioFormatsImage", function(x, v, ...) {
    ext(x) <- .shift_ext(ext(x), v)
    origin(x) <- origin(x) + v
    x
})

#' @rdname translateImg
#' @export
setMethod("translateImg", "ExtImage", function(x, v, ...) {
    ext(x) <- .shift_ext(ext(x), v)
    x
})

# Scale------------------
.scale_ext <- function(x, factor, bbox_all = ext(x), ...) {
    ext(x) <- .transform_bbox(ext(x), list(name="scale", factor=factor),
                              bbox_all = bbox_all)
    x
}
#' Scale image
#'
#' This function scales the image about its center. After scaling, the center
#' of the image is not shifted.
#'
#' @inheritParams transposeImg
#' @param factor Numeric, scaling factor.
#' @return A \code{*Image} object of the same class that has been scaled. Behind
#' the scene, it's only the extent that has been changed and the images are not
#' changed. The center of the image is unchanged.
#' @aliases scaleImg
#' @name scaleImg
#' @concept Image affine transformation
#' @family image methods
#' @export
NULL

#' @rdname scaleImg
#' @export
setMethod("scaleImg", "AlignedSpatialImage",
          function(x, factor, ...) .scale_ext(x, factor))

# Affine ------------

#' Affine transformation of images
#'
#' This function performs affine transformation on images, with any matrix and
#' translation vector.
#'
#' @inheritParams transposeImg
#' @param M A 2x2 numeric matrix for the linear transformation in the xy plane.
#' @param v A numeric vector of length 2 for translation in the xy plane.
#' @return \code{SpatRasterImage} will be converted to \code{ExtImage}. Otherwise
#' \code{*Image} object of the same class. For \code{BioFormatsImage}, the
#' transformation info is stored and will be applied when the image is loaded
#' into memory as \code{ExtImage}.
#' @name affineImg
#' @aliases affineImg
#' @concept Image affine transformation
#' @family image methods
#' @export
NULL

#' @rdname affineImg
#' @export
setMethod("affineImg", "SpatRasterImage", # Deal with rotating SpatRaster?
          function(x, M, v, maxcell = 1e7, ...) {
              # Not sure what exactly to do. I think convert to ExtImage as well.
              # Covert to BioFormatsImage if it's not in memory
              x <- toExtImage(x, maxcell)
              affineImg(x, M, v)
          })

#' @rdname affineImg
#' @export
setMethod("affineImg", "BioFormatsImage",
          function(x, M, v, ...) {
              .combine_transforms(x, list(M=M, v=v))
          })

#' @rdname affineImg
#' @export
setMethod("affineImg", "ExtImage",
          function(x, M, v, ...) {
              .affine_ebi(x, ext(x), M, v)
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
#' @concept Image methods
#' @export
#' @family image methods
NULL

#' @rdname cropImg
#' @export
setMethod("cropImg", "SpatRasterImage", function(x, bbox, filename = "") {
    x <- terra::crop(x, bbox, snap = "out", filename = filename)
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
setMethod("cropImg", "ExtImage", function(x, bbox) {
    # Convert bbox to pixel range based on ext(x)
    bbox_old <- ext(x)
    origin <- bbox_old[c("xmin", "ymin")]
    bbox <- .shift_ext(bbox, -origin)
    dim_old <- dim(x)
    sfx <- dim_old[1]/(bbox_old["xmax"] - bbox_old["xmin"])
    sfy <- dim_old[2]/(bbox_old["ymax"] - bbox_old["ymin"])
    bbox_use <- bbox
    bbox_use[c("xmin", "xmax")] <- bbox_use[c("xmin", "xmax")] * sfx
    bbox_use[c("ymin", "ymax")] <- bbox_use[c("ymin", "ymax")] * sfy
    min_nms <- c("xmin", "ymin")
    max_nms <- c("xmax", "ymax")

    subset_use <- .get_subset(bbox_use, dim_old[2])
    bbox_use[min_nms] <- floor(bbox_use[min_nms])
    bbox_use[max_nms] <- ceiling(bbox_use[max_nms])

    if (length(dim_old) == 3L)
        x <- x[subset_use[[1]], subset_use[[2]],]
    else
        x <- x[subset_use[[1]], subset_use[[2]]]

    bbox_new <- bbox_use
    bbox_new[c("xmin", "xmax")] <- bbox_new[c("xmin", "xmax")] / sfx
    bbox_new[c("ymin", "ymax")] <- bbox_new[c("ymin", "ymax")] / sfy
    ext(x) <- .shift_ext(bbox_new, origin)
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
