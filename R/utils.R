#' Get all unique sample IDs
#'
#' The title is self-explanatory.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @return A character vector of all unique entries of the \code{sample_id}
#' column in \code{colData(x)}.
#' @export
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @concept Utilities
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sampleIDs(sfe)
sampleIDs <- function(sfe) unique(colData(sfe)$sample_id)

#' Change sample IDs
#'
#' Change sample IDs in all fields of the SFE object where sample IDs are
#' present, not just the colData.
#'
#' @inheritParams sampleIDs
#' @param replacement A named character vector whose names are the existing
#'   sample IDs to be changed and whose values are the corresponding
#'   replacements.
#' @return An SFE object.
#' @concept Utilities
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
#' sampleIDs(sfe)
changeSampleIDs <- function(sfe, replacement) {
    for (i in seq_along(replacement)) {
        original <- names(replacement)[i]
        colData(sfe)$sample_id[colData(sfe)$sample_id == original] <-
            replacement[i]
        gs_names <- names(int_metadata(sfe)$spatialGraphs)
        names(int_metadata(sfe)$spatialGraphs)[gs_names == original] <-
            replacement[i]
        if (length(int_metadata(sfe)$annotGeometries)) {
            for (n in names(int_metadata(sfe)$annotGeometries)) {
                ag <- int_metadata(sfe)$annotGeometries[[n]]
                ind <- ag$sample_id == original
                int_metadata(sfe)$annotGeometries[[n]]$sample_id[ind] <- replacement[i]
            }
        }
        if (length(rowGeometries(sfe))) {
            nms <- rowGeometryNames(sfe)
            nms <- gsub(paste0(original, "$"), replacement[i], nms)
            rowGeometryNames(sfe) <- nms
            # Edge case: what if one sample_id includes another one?
            # e.g. sample01_x and x
        }
        if (nrow(imgData(sfe))) {
            imgData(sfe)$sample_id[imgData(sfe)$sample_id == original] <-
                replacement[i]
        }
        # Check spatial results
        # rowData
        rd_ind <- grepl(paste0("_", original), names(rowData(sfe)))
        if (any(rd_ind)) {
            names(rowData(sfe)) <- gsub(paste0("_", original), paste0("_", replacement[i]),
                                        names(rowData(sfe)))
        }
        # featureData
        if (!is.null(colFeatureData(sfe))) {
            nms <- names(colFeatureData(sfe))
            nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
            names(colFeatureData(sfe)) <- nms
        }
        # colGeometries
        for (n in colGeometryNames(sfe)) {
            if (!is.null(geometryFeatureData(sfe, n, 2L))) {
                nms <- names(geometryFeatureData(sfe, n, 2L))
                nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
                names(geometryFeatureData(sfe, n, 2L)) <- nms
            }
        }
        # annotGeometries
        for (n in annotGeometryNames(sfe)) {
            if (!is.null(geometryFeatureData(sfe, n, 3L))) {
                nms <- names(geometryFeatureData(sfe, n, 3L))
                nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
                names(geometryFeatureData(sfe, n, 3L)) <- nms
            }
        }
        # reducedDims
        for (n in reducedDimNames(sfe)) {
            if (!is.null(reducedDimFeatureData(sfe, n))) {
                nms <- names(reducedDimFeatureData(sfe, n))
                nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
                names(reducedDimFeatureData(sfe, n)) <- nms
            }
        }
    }
    sfe
}

.translate_value <- function(x, translate, value, sample_id = NULL) {
    if (translate && !is.null(int_metadata(x)$orig_bbox)) {
        if (anyNA(value$sample_id) || is.null(value$sample_id)) {
            if (nrow(value) == ncol(x))
                value$sample_id <- colData(x)$sample_id
            else if (nrow(value) == nrow(x))
                value$sample_id <- .check_sample_id(x, sample_id)
        }
        orig_bbox <- int_metadata(x)$orig_bbox
        # Don't translate if already translated
        curr_bbox <- st_bbox(value)
        samples <- unique(value$sample_id) %||% sample_id
        if (length(samples) > 1L) {
            value$ID_ <- seq_len(nrow(value)) # Unlikely name
            df <- value[,c("ID_", "sample_id", "geometry")]
            df_split <- split(df, value$sample_id)
            df_split <- lapply(samples, function(s) {
                out <- df_split[[s]]
                og <- out$geometry - orig_bbox[c("xmin", "ymin"), s]
                bb <- st_bbox(og)
                if ((bb["xmin"] < 0 || bb["ymin"] < 0) || anyNA(bb)) return(out)
                else
                    out$geometry <- og
                out
            })
            df <- do.call(rbind, df_split)
            df <- df[match(value$ID_, df$ID_),]
            value$geometry <- df$geometry
            value$ID_ <- NULL
        } else {
            if (curr_bbox["xmin"] - orig_bbox["xmin", samples] < 0 ||
                curr_bbox["ymin"] - orig_bbox["ymin", samples] < 0) {
                return(value)
            }
            value$geometry <- value$geometry - orig_bbox[c("xmin", "ymin"), samples]
        }
    }
    value
}

.path_valid2 <- function(x) {
    all(c(length(x) == 1, is.character(x), file.exists(x)))
}

#' Find center of bounding box
#'
#' Get x-y coordinates of the center of any bounding box
#'
#' @param bbox A numeric vector of length 4 with names xmin, xmax, ymin, ymax,
#' in any order.
#' @return A numeric vector of length 2.
#' @concept Utilities
#' @examples
#' bbox <- c(xmin = 0, xmax = 100, ymin = 0, ymax = 80)
#' bbox_center(bbox)
#' @export
bbox_center <- function(bbox) {
    c(mean(bbox[c("xmin", "xmax")]),
      mean(bbox[c("ymin", "ymax")]))
}

#' Get physical size of pixels
#'
#' This function gets physical size of pixels in each resolution of a OME-TIFF
#' pyramid in \code{\link{BioFormatsImage}}.
#'
#' @param file Path to an OME-TIFF file.
#' @param resolution Which resolution to query; 1 means the highest resolution.
#' The pixels will be larger for the lower resolutions.
#' @return Numeric vector of length 2 of pixel size in x and y. Usually they're
#' the same.
#' @export
#' @concept Utilities
#' @examples
#' library(SFEData)
#' fp <- tempdir()
#' dir_use <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
#' # RBioFormats null pointer error
#' try(getPixelSize(file.path(dir_use, "morphology_focus.ome.tif")))
#' getPixelSize(file.path(dir_use, "morphology_focus.ome.tif"))
#' unlink(dir_use, recursive = TRUE)
getPixelSize <- function(file, resolution = 1L) {
    check_installed(c("xml2", "RBioFormats"))
    xml_meta <- RBioFormats::read.omexml(file) |>
        xml2::read_xml() |> xml2::as_list()
    psx <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeX") |> as.numeric()
    psy <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeY") |> as.numeric()
    if (resolution == 1L) return(c(psx, psy))
    else {
        m <- RBioFormats::read.metadata(file)
        coreMetadata <- RBioFormats::coreMetadata
        meta <- coreMetadata(m, series = resolution)
        meta1 <- coreMetadata(m, series = 1L)
        sizeX_full <- meta1$sizeX
        sizeY_full <- meta1$sizeY
        fct_x <- sizeX_full/meta$sizeX
        fct_y <- sizeY_full/meta$sizeY
        fct_round <- round(fct_x) # Should be the same for x and y
        fctx2 <- fct_x/fct_round
        fcty2 <- fct_y/fct_round

        sfx2 <- meta$sizeX*fctx2/sizeX_full # Multiply to get from full res pixel to low res pixel
        sfy2 <- meta$sizeY*fcty2/sizeY_full
        psx_out <- psx/sfx2
        psy_out <- psy/sfy2
        return(c(psx_out, psy_out))
    }
}

#' Aggregate bounding boxes
#'
#' To find the bounding box of multiple bounding boxes.
#'
#' @param bboxes Either a matrix with 4 rows whose columns are the different
#' bounding boxes, with row names "xmin", "xmax", "ymin", and "ymax" in any order,
#' or a list of bounding boxes which are named numeric vectors.
#' @return A named numeric vector for the total bounding box.
#' @export
#' @concept Utilities
#' @examples
#' bboxes <- list(c(xmin = 5, xmax = 10, ymin = 2, ymax = 20),
#' c(xmin = 8, xmax = 18, ymin = 0, ymax = 15))
#' bbox_all <- aggBboxes(bboxes)
aggBboxes <- function(bboxes) {
    if (is.list(bboxes)) {
        bboxes <- lapply(bboxes, function(x) x[c("xmin", "ymin", "xmax", "ymax")])
        bboxes <- do.call(rbind, bboxes)
    }
    c(
        xmin = min(bboxes[, "xmin"], na.rm = TRUE),
        ymin = min(bboxes[, "ymin"], na.rm = TRUE),
        xmax = max(bboxes[, "xmax"], na.rm = TRUE),
        ymax = max(bboxes[, "ymax"], na.rm = TRUE)
    )
}

#' Show all image_ids in the SFE object
#'
#' The title is self-explanatory. Some functions require \code{image_id} to get
#' or set images.
#'
#' @inheritParams sampleIDs
#' @return A character vector of \code{image_ids}.
#' @export
#' @concept Utilities
#' @examples
#' fp <- system.file(file.path("extdata", "sample01"),
#' package = "SpatialFeatureExperiment")
#' sfe <- read10xVisiumSFE(fp, type = "sparse")
#' imageIDs(sfe)
imageIDs <- function(sfe) imgData(sfe)$image_id

#' Check if Parquet GDAL driver is available
#'
#' The GeoParquet files for geometries are typically written and read with the
#' \code{sfarrow} package, but to add only a select few genes to the SFE object
#' say for visualization purposes, the Parquet GDAL driver is required in order
#' to use GDAL's SQL to query the GeoParquet file to only load the few genes
#' requested. The transcript spots from a large dataset can take up a lot of
#' memory if all loaded.
#'
#' The Parquet driver has been supported since GDAL 3.5.0. The \code{arrow} C++
#' library must be installed in order to make the Parquet driver available. When
#' arrow is installed, newer versions of GDAL installed from Homebrew (Mac)
#' should have the Parquet driver. For Linux, the binary from \code{apt-get}'s
#' default repo is 3.4.1 (as of April 2024). To use the Parquet driver, GDAL may
#' need to be installed from source. See script from the \href{https://github.com/rocker-org/rocker-versioned2/blob/master/scripts/experimental/install_dev_osgeo.sh}{geospatial rocker}.
#' A Voyager docker container with the Parquet driver will soon be provided.
#'
#' @return Logical, indicating whether the Parquet driver is present.
#' @export
#' @concept Utilities
#' @examples
#' gdalParquetAvailable()
#'
gdalParquetAvailable <- function() {
    "Parquet" %in% rownames(sf::st_drivers())
}

.get_XOA_version <- function(data_dir) {
    experiment <- fromJSON(file = file.path(data_dir, "experiment.xenium"),
                           simplify = TRUE)
    xoa_version <- experiment$analysis_sw_version
    major_version <- substr(xoa_version, 8, 8) |> as.integer()
    minor_version <- substr(xoa_version, 10, 10) |> as.integer()
    instrument_version <- experiment$instrument_sw_version
    c(xoa = xoa_version, major = major_version, minor = minor_version,
      instrument = instrument_version)
}

.no_raw_bytes <- function(data_dir) {
    c(xoa_version, major_version, minor_version, instrument_version) %<-%
        .get_XOA_version(data_dir)
    (major_version == 1L && minor_version > 4L) || major_version == 2L
}

.get_tech_tx_fields <- function(tech, data_dir) {
    spatialCoordsNames <- switch(
        tech,
        Vizgen = c("global_x", "global_y", "global_z"),
        Xenium = c("x_location", "y_location", "z_location"),
        CosMX = c("x_global_px", "y_global_px", "z")
    )
    gene_col <- switch(
        tech,
        CosMX = "target",
        Xenium = "feature_name",
        Vizgen = "gene"
    )
    cell_col <- switch(
        tech,
        Vizgen = "barcode_id",
        Xenium = "cell_id",
        CosMX = "cell_ID"
    )
    fn <- switch(
        tech,
        Vizgen = .check_vizgen_fns(data_dir, "detected_transcripts.csv"),
        CosMX = grep("tx_file.csv",
                     list.files(data_dir, pattern = "\\.csv$", full.names = TRUE),
                     value = TRUE),
        Xenium = .check_xenium_fns(data_dir, "transcripts", .no_raw_bytes(data_dir))
    )
    list(spatialCoordsNames = spatialCoordsNames,
         gene_col = gene_col,
         cell_col = cell_col,
         fn = fn)
}

.id2symbol <- function(x, ids, swap_rownames) {
    if (!is.null(swap_rownames) && swap_rownames %in% names(rowData(x))) {
        inds <- match(ids, rownames(x))
        features <- rowData(x)[[swap_rownames]][inds]
    } else features <- ids
    features
}
