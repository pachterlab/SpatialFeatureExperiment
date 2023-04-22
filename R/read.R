#' Read 10X Visium data as SpatialFeatureExperiment
#'
#' Read Space Ranger output as a SpatialFeatureExperiment object, where spots
#' are represented with polygons in the colGeometry called "spotPoly". Other
#' geometries can be added later after the dataset is read. If \code{data =
#' "filtered"}, then spatial neighborhood graphs of the spots are also computed
#' and stored in the colGraph called "visium" in all samples for downstream
#' spatial analyses.
#'
#' @inheritParams SpatialExperiment::read10xVisium
#' @inheritParams findVisiumGraph
#' @inheritParams SpatialFeatureExperiment
#' @param type Either "HDF5", and the matrix will be represented as
#'   \code{TENxMatrix}, or "sparse", and the matrix will be read as a
#'   \code{dgCMatrix}.
#' @param dirs Directory for each sample that contains the \code{spatial} and
#'   \code{raw/filtered_featues_bc_matrix} directories. By default, the
#'   \code{outs} directory under the directory specified in the \code{samples}
#'   argument, as in Space Ranger output. Change the \code{dirs} argument if you
#'   have moved or renamed the output directory.
#' @param unit Whether to use pixels in full resolution image or microns as the
#'   unit. If using microns, then spacing between spots in pixels will be used
#'   to convert the coordinates into microns, as the spacing is known to be 100
#'   microns. This is used to plot scale bar.
#' @param load Not used, kept for backward compatibility.
#' @note The \code{as(<dgTMatrix>, "dgCMatrix") is deprecated} warning comes
#'   from the \code{DropletUtils} package which is used by
#'   \code{SpatialExperiment} to read 10X outputs. This will be fixed when
#'   \code{SpatialExperiment} switches to TENxIO.
#' @importFrom SpatialExperiment read10xVisium
#' @importFrom rjson fromJSON
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom utils read.csv
#' @note It is assumed that the images have not been cropped. Otherwise the
#' images might not align with the spots.
#' @return A SpatialFeatureExperiment object. The images might need to be
#'   manually transposed and/or mirrored to match the spots in this version of
#'   this package.
#' @export
#' @examples
#' dir <- system.file("extdata", package = "SpatialFeatureExperiment")
#'
#' sample_ids <- c("sample01", "sample02")
#' samples <- file.path(dir, sample_ids)
#'
#' list.files(samples[1])
#' list.files(file.path(samples[1], "spatial"))
#' (sfe <- read10xVisiumSFE(samples, sample_id = sample_ids,
#'     type = "sparse", data = "filtered",
#'     load = FALSE
#' ))
read10xVisiumSFE <- function(samples = "",
                             dirs = file.path(samples, "outs"),
                             sample_id = paste0(
                                 "sample",
                                 sprintf(
                                     "%02d",
                                     seq_along(samples)
                                 )
                             ),
                             type = c("HDF5", "sparse"),
                             data = c("filtered", "raw"),
                             images = c("lowres", "hires"),
                             unit = c("full_res_image_pixel", "micron"),
                             style = "W", zero.policy = NULL,
                             BPPARAM = SerialParam(), load = FALSE) {
    type <- match.arg(type)
    data <- match.arg(data)
    unit <- match.arg(unit)
    images <- match.arg(images, several.ok = TRUE)
    img_fns <- c(
        lowres="tissue_lowres_image.png",
        hires="tissue_hires_image.png")
    img_fns <- img_fns[images]
    # Read one sample at a time, in order to get spot diameter one sample at a time
    sfes <- lapply(seq_along(samples), function(i) {
        o <- read10xVisium(dirs[i], sample_id[i], type, data, images, load = FALSE)
        imgData(o) <- NULL

        scalefactors <- fromJSON(file = file.path(
            dirs[i], "spatial",
            "scalefactors_json.json"
        ))
        o <- .spe_to_sfe(o,
            colGeometries = NULL, rowGeometries = NULL,
            annotGeometries = NULL, spatialCoordsNames = NULL,
            annotGeometryType = NULL, spatialGraphs = NULL,
            spotDiameter = scalefactors[["spot_diameter_fullres"]],
            unit = unit, BPPARAM = BPPARAM
        )
        # Add spatial enrichment if present
        fn <- file.path(dirs[i], "spatial", "spatial_enrichment.csv")
        if (file.exists(fn)) {
            enrichment <- read.csv(fn)
            row_inds <- match(rownames(o), enrichment$Feature.ID)
            # Let me not worry about different samples having different genes for now
            if (i == 1L) {
                rowData(o) <- cbind(rowData(o),
                                    Feature.Type = enrichment[row_inds, "Feature.Type"])
            }
            cols_use <- c("I", "P.value", "Adjusted.p.value",
                          "Feature.Counts.in.Spots.Under.Tissue",
                          "Median.Normalized.Average.Counts",
                          "Barcodes.Detected.per.Feature")
            enrichment2 <- enrichment[row_inds, cols_use]
            names(enrichment2) <- paste(names(enrichment2), sample_id[i],
                                        sep = "_")
            rowData(o) <- cbind(rowData(o), enrichment2)
        }
        # Add barcode fluorescence intensity if present
        fn2 <- file.path(dirs[i], "spatial", "barcode_fluorescence_intensity.csv")
        if (file.exists(fn2)) {
            fluo <- read.csv(fn2)
            row_inds <- match(colnames(o), fluo$barcode)
            fluo$barcode <- NULL
            fluo$in_tissue <- NULL
            colData(o) <- cbind(colData(o), fluo[row_inds,])
        }

        names_use <- paste("tissue", images, "scalef", sep = "_")
        scale_imgs <- unlist(scalefactors[names_use])
        # Convert to microns and set extent for image
        if (unit == "micron") {
            scale_fct <- .pixel2micron(o)
            cg <- spotPoly(o)
            cg$geometry <- cg$geometry * scale_fct
            spotPoly(o) <- cg
            # Scale factors for images
            scale_imgs <- scale_imgs / scale_fct
        } else {
            scale_imgs <- scalefactors[paste("tissue", images, "scalef",
                                             sep = "_")]
        }
        # Set up ImgData
        img_fns2 <- file.path(dirs[i], "spatial", img_fns)
        img_dfs <- lapply(seq_along(img_fns), function(j) {
            .get_imgData(img_fns2[j], sample_id = sample_id[i],
                         image_id = names(img_fns)[j],
                         extent = NULL, scale_fct = scale_imgs[[j]],
                         flip = TRUE)
        })
        img_df <- do.call(rbind, img_dfs)
        imgData(o) <- img_df
        o
    })
    out <- do.call(cbind, sfes)
    if (data == "filtered") {
        colGraphs(out, sample_id = "all", name = "visium") <-
            findVisiumGraph(out,
                sample_id = "all", style = style,
                zero.policy = zero.policy
            )
    }
    out
}

#' @importFrom sf st_nearest_feature st_distance
.pixel2micron <- function(sfe) {
    min_row <- min(sfe$array_row)
    min_col <- min(sfe$array_col)
    inds_sub <- sfe$array_row <= min_row + 8 & sfe$array_col <= min_col + 8
    coords_sub <- df2sf(spatialCoords(sfe)[inds_sub,], spatialCoordsNames(sfe))
    inds <- st_nearest_feature(coords_sub)
    dists <- vapply(seq_along(inds), function(i) {
        st_distance(coords_sub[i,], coords_sub[inds[i],])[1,1]
    }, FUN.VALUE = numeric(1))
    dists <- mean(dists) # Full res pixels per 100 microns
    100/dists
}

.h52poly_fov <- function(fn, z) {
    l <- rhdf5::h5dump(fn)[[1]]
    cell_ids <- names(l)
    z_name <- paste0("zIndex_", z)
    geometries <- lapply(l, function(m) {
        sf::st_polygon(list(t(m[[z_name]]$p_0$coordinates[,,1])))
    })
    df <- data.frame(geometry = sf::st_sfc(geometries),
                     ID = cell_ids,
                     ZIndex = z)
    sf::st_sf(df)
}

.filter_polygons <- function(polys, min_area) {
    if (st_geometry_type(polys, by_geometry = FALSE) == "MULTIPOLYGON") {
        polys_sep <- lapply(st_geometry(polys), function(x) {
            st_cast(st_sfc(x), "POLYGON")
        })
        areas <- lapply(polys_sep, st_area)
        which_keep <- lapply(areas, function(x) which(x > min_area))
        multi_inds <- which(lengths(which_keep) > 1L)
        if (length(multi_inds)) {
            warning("There are ", length(multi_inds), " cells with multiple",
                    " pieces in cell segmentation larger than min_area,",
                    " whose indices are: ",
                    paste(multi_inds, collapse = ", "),
                    ". The largest piece is kept.")
            which_keep[multi_inds] <- lapply(areas[multi_inds], which.max)
        }
        inds <- lengths(which_keep) > 0L
        polys <- polys[inds,]
        which_keep <- unlist(which_keep[inds])
        geo <- st_geometry(polys)
        new_geo <- lapply(seq_along(which_keep), function(i) {
            geo[[i]] <- st_cast(geo[i], "POLYGON")[[which_keep[i]]]
        }) |> st_sfc()
        st_geometry(polys) <- new_geo
    } else {
        inds <- st_area(st_geometry(polys)) > min_area
        polys <- polys[inds,]
    }
    polys
}

.if_flip_img <- function(fn, max_flip) {
    max_flip <- toupper(max_flip)
    unit <- gsub("^[0-9.]+\\s?", "", max_flip)
    if (!unit %in% c("MB", "GB"))
        stop("max_flip must be in either MB or GB.")
    max_flip <- as.numeric(gsub("\\s?[A-Z.]+$", "", max_flip))
    max_flip <- switch (unit,
        MB = max_flip * 1024^2,
        GB = max_flip * 1024^3
    )
    size <- file.info(fn)[["size"]] # NA if file doesn't exist
    size < max_flip
}

#' Read Vizgen MERFISH output as SpatialFeatureExperiment
#'
#' This function reads the standard Vizgen MERFISH output into an SFE object.
#' The coordinates are in microns. Cell centroids are read into
#' \code{\link{colGeometry}} "centroids", and cell segmentations are read into
#' \code{colGeometry} "cellSeg". The image(s) (polyT and DAPI) are also read as
#' \code{\link{SpatRaster}} objects so they are not loaded into memory unless
#' necessary. Because the image's origin is the top left while the geometry's
#' origin is bottom left, either the image or the geometry needs to be flipped.
#' Because the image accompanying MERFISH datasets are usually very large, the
#' coordinates will be flipped so the flipping operation won't load the entire
#' image into memory.
#'
#' @inheritParams SpatialFeatureExperiment
#' @param data_dir Top level directory of Vizgen output, which contains
#'   directories \code{cell_boundaries} and \code{images}, and files
#'   \code{cell_by_gene.csv}, \code{cell_metadata.csv}, and
#'   \code{detected_transcripts.csv}.
#' @param z Index of z plane to read.
#' @param use_cellpose Logical, whether to use Cellpose parquet files if
#'   present.
#' @param max_flip Maximum size of the image allowed to flip the image. Because
#' the image will be loaded into memory to be flipped. If the image is larger
#' than this size then the coordinates will be flipped instead.
#' @param flip To flip the image, geometry coordinates, or none. Because the
#'   image has the origin at the top left while the geometry has origin at the
#'   bottom left, one of them needs to be flipped for them to match. If one of
#'   them is already flipped, then use "none". The image will not be flipped if
#'   it's GeoTIFF.
#' @param image Which image(s) to load, can be "DAPI", "PolyT", or both.
#' @param min_area Minimum cell area in square microns. Anything smaller will be
#'   considered artifact or debris and removed.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying parallel
#'   processing backend and number of threads to use to load cell segmentation
#'   from HDF5 files from different fields of view (FOVs) with multiple cores. A
#'   progress bar can be configured in the \code{\link{BiocParallelParam}}
#'   object. When there are numerous FOVs, reading in the geometries can be time
#'   consuming, so we recommend using a server and larger number of threads.
#'   This argument is not used if \code{use_cellpose = TRUE} and the parquet
#'   file is present.
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
#' @importFrom sf st_area st_geometry<-
#' @importFrom terra rast ext vect
#' @importFrom BiocParallel bpmapply
#' @importFrom rlang check_installed
#' @examples
#' dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
#' sfe <- readVizgen(dir_use, z = 0L, use_cellpose = TRUE, image = "PolyT",
#' flip = "geometry")
readVizgen <- function(data_dir, z = 3L, use_cellpose = TRUE,
                       sample_id = "sample01", min_area = 15,
                       image = c("DAPI", "PolyT"),
                       flip = c("geometry", "image", "none"),
                       max_flip = "50 MB",
                       BPPARAM = SerialParam()) {
    data_dir <- normalizePath(data_dir, mustWork = TRUE)
    flip <- match.arg(flip)
    image <- match.arg(image, several.ok = TRUE)
    rlang::check_installed("vroom")
    if (z < 0 || z > 6) {
        stop("z must be beween 0 and 6 (inclusive).")
    }
    img_fn <- file.path(data_dir, "images", paste0("mosaic_", image, "_z", z, ".tif"))
    names(img_fn) <- image
    do_flip <- .if_flip_img(img_fn, max_flip)
    if_exists <- file.exists(img_fn)
    if (!all(if_exists)) {
        warning("The image file(s) for ", image[!if_exists],
                " don't exist, or have non-standard file name(s).")
        img_fn <- img_fn[if_exists]
    }
    if (!length(img_fn)) flip <- "none"
    else if (!any(do_flip) && flip == "image") flip <- "geometry"
    parq_files <- list.files(data_dir, "*.parquet")
    use_cellpose <- use_cellpose & length(parq_files)
    if (use_cellpose) {
        rlang::check_installed("sfarrow")
        fn <- file.path(data_dir, "cellpose_micron_space.parquet")
        polys <- sfarrow::st_read_parquet(fn)
        polys <- polys[polys$ZIndex == z,]
        polys <- .filter_polygons(polys, min_area)
        st_geometry(polys) <- "geometry"
        polys$ID <- polys$EntityID
        polys <- polys[,c("ID", "ZIndex", "Type", "ZLevel", "geometry")]
    } else {
        rlang::check_installed("rhdf5")
        rlang::check_installed("dplyr")
        fns <- list.files(file.path(data_dir, "cell_boundaries"),
                          "*.hdf5", full.names = TRUE)
        polys <- bpmapply(.h52poly_fov, fn = fns, SIMPLIFY = FALSE,
                          BPPARAM = BPPARAM,
                          MoreArgs = list(z = z))
        # dplyr is much more efficient than base R rbind
        polys <- if (length(polys) == 1L) polys[[1]] else do.call(dplyr::bind_rows, polys)
    }
    if (flip == "geometry") {
        # Flip the coordinates
        mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
        st_geometry(polys) <- st_geometry(polys) * mat_flip
    }

    mat_fn <- list.files(data_dir, "cell_by_gene.csv", full.names = TRUE)
    suppressMessages(mat <- vroom::vroom(mat_fn, col_types = vroom::cols(...1 = "c")))
    mat <- mat[match(polys$ID, mat$...1),]
    meta_fn <- list.files(data_dir, "cell_metadata.csv", full.names = TRUE)
    suppressMessages(metadata <- vroom::vroom(meta_fn, col_types = vroom::cols(...1 = "c")))
    metadata <- metadata[match(polys$ID, metadata$...1),]
    if (flip == "geometry") {
        metadata$center_y <- -metadata$center_y
    }

    m <- as.matrix(mat[,-1])
    m <- as(m, "CsparseMatrix")
    rownames(m) <- mat$...1
    m <- Matrix::t(m)

    manifest <- fromJSON(file = file.path(data_dir, "images", "manifest.json"))
    extent <- setNames(manifest$bbox_microns, c("xmin", "ymin", "xmax", "ymax"))
    if (flip == "geometry") {
        extent[c("ymin", "ymax")] <- -extent[c("ymax", "ymin")]
    }
    # Set up ImgData
    img_dfs <- lapply(names(img_fn), function(n) {
        .get_imgData(img_fn[n], sample_id = sample_id, image_id = n,
                     extent = extent,
                     flip = (flip == "image"))
    })
    img_df <- do.call(rbind, img_dfs)
    # Takes a while to make the POINT geometry for the centroids, not too bad
    sfe <- SpatialFeatureExperiment(assays = list(counts = m),
                                    colData = metadata[,-1],
                                    spatialCoordsNames = c("center_x", "center_y"),
                                    unit = "micron")
    rownames(polys) <- polys$ID
    polys$ID <- NULL
    cellSeg(sfe) <- polys
    imgData(sfe) <- img_df
    sfe
}
