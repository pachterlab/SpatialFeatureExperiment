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
#' @concept Read data into SFE
#' @importFrom DropletUtils read10xCounts
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
                             style = "W", zero.policy = NULL, load = FALSE) {
    type <- match.arg(type)
    data <- match.arg(data)
    unit <- match.arg(unit)
    images <- match.arg(images, several.ok = TRUE)
    img_fns <- c(
        lowres="tissue_lowres_image.png",
        hires="tissue_hires_image.png")
    img_fns <- img_fns[images]
    
    # Read one sample at a time, in order to get spot diameter one sample at a time
    # TODO: need support for VisiumHD ----
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
                         unit = unit
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
        # Create Visium graph for filtered data
        if (data == "filtered") {
            colGraph(o, "visium") <- findVisiumGraph(
                o, sample_id = "all", style = style,
                zero.policy = zero.policy
            )
        }
        o
    })
    out <- do.call(cbind, sfes)
    out
}

#' @importFrom sf st_nearest_feature st_distance
#' @importFrom stats median
.pixel2micron <- function(sfe) {
    # Use center spots rather than corner, to be more robust for filtered data
    mid_row <- median(sfe$array_row)
    mid_col <- median(sfe$array_col)
    inds_sub <- abs(sfe$array_row - mid_row) <= 2 & abs(sfe$array_col - mid_col) <= 2
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
        m[[z_name]]$p_0$coordinates[, , 1]
    })
    # NOTE: some coordinates can be empty, ie NULL
    # get boolean indices of cells, ie TRUE for non-empty
    inds <- lengths(geometries) > 0
    # remove empty elements
    geometries <- geometries[inds]
    geometries <- lapply(geometries, function(m) st_polygon(list(t(m))))

    # keep non-emplty elements
    df <- st_sf(geometry = sf::st_sfc(geometries),
                ID = cell_ids[which(inds)],
                ZIndex = z)
    df
}

#' @importFrom sf st_is_empty
#' @importFrom BiocParallel bplapply
#' @importFrom utils head
.filter_polygons <- function(polys, min_area, BPPARAM = SerialParam()) {
    # Sanity check on nested polygon lists
    test.segs <- vapply(st_geometry(polys), length, FUN.VALUE = integer(1))
    if (any(test.segs > 1)) {
        segs.art.index <- which(test.segs > 1)
        warning("Sanity checks on cell segmentation polygons:", "\n",
                ">>> ..found ", length(segs.art.index),
                " cells with (nested) polygon lists", "\n",
                ">>> ..applying filtering")
    }
    # remove empty elements
    polys.orig <- polys
    polys <- polys[!st_is_empty(polys),]
    empty.inds <- which(!polys.orig$ID %in% polys$ID)
    if (length(empty.inds)) { message(">>> ..removing ",
                                      length(empty.inds), " empty polygons") }
    if (st_geometry_type(polys, by_geometry = FALSE) == "MULTIPOLYGON") {
        polys_sep <- lapply(st_geometry(polys), function(x) {
            st_cast(st_sfc(x), "POLYGON")
        })
        areas <- lapply(polys_sep, st_area)

        if (!is.null(min_area)) {
            which_keep <- lapply(areas, function(x) which(x > min_area))
            multi_inds <- which(lengths(which_keep) > 1L)
            if (length(multi_inds)) {
                warning("There are ", length(multi_inds), " cells with multiple",
                        " pieces in cell segmentation larger than min_area,",
                        " whose first 10 indices are: ",
                        paste(multi_inds |> head(10), # necessary to print all?
                              collapse = ", "),
                        ". The largest piece is kept.")
                which_keep[multi_inds] <- lapply(areas[multi_inds], which.max)
            }
            inds <- lengths(which_keep) > 0L
            polys <- polys[inds,]
            # using parallelization, else can take a while when `which_keep` length is towards 100K
            which_keep <- unlist(which_keep[inds])
        } else if (is.null(min_area)) {
            # use only maximal area, ie the largest polygon
            warning(">>> ..keeping polygons with the largest area only")
            which_keep <- lapply(areas, function(x) which.max(x))
        }
        geo <- st_geometry(polys)
        new_geo <- # Should only iterate over those with multiple pieces
            bplapply(seq_along(which_keep), function(i) {
                geo[[i]] <- st_cast(geo[i], "POLYGON")[[which_keep[i]]] |>
                    unique() |> # remove any duplicates
                    st_polygon()
            }, BPPARAM = BPPARAM) |> st_sfc()
        st_geometry(polys) <- new_geo
    } else {
        inds <- st_area(st_geometry(polys)) > min_area
        if (any(inds)) {
            message("Removing ", sum(!inds), " cells with area less than ", min_area)
        }
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

.check_vizgen_fns <- function(data_dir, keyword) {
    fn <-
        list.files(data_dir,
                   pattern = keyword,
                   full.names = TRUE)
    if (length(fn) == 0) {
        fn <- list.files(data_dir,
                         pattern = keyword,
                         full.names = TRUE,
                         recursive = TRUE)
    }
    if (length(fn) > 1) {
        stop("There are > 1 `", keyword, "` files",
             "\n", "make sure only 1 file is read")
    } else if (!length(fn)) {
        stop("No `", keyword, "` file is available")
    }
    fn
}

# sanity on geometries to remove any self-intersection
#' @importFrom sf st_buffer st_is_valid
.check_st_valid <- function(sf_df = NULL) {
    # sf_df is a single sf data frame, not a list of sf data frames
    invalid_inds <- which(!st_is_valid(sf_df))
    if (length(invalid_inds)) {
        geoms_new <- st_buffer(st_geometry(sf_df)[invalid_inds], dist = 0)
        # [.sfc<- is slow for st_GEOMETRY when buffer created MULTIPOLYGONs
        # due to a vapply somewhere to recompute bbox in an R loop
        # But not too bad, I don't think worse than st_cast
        st_geometry(sf_df)[invalid_inds] <- geoms_new
    }
    return(sf_df)
}

#' Read Vizgen MERFISH output as SpatialFeatureExperiment
#'
#' This function reads the standard Vizgen MERFISH output into an SFE object.
#' The coordinates are in microns. Cell centroids are read into
#' \code{\link{colGeometry}} "centroids", and cell segmentations are read into
#' \code{colGeometry} "cellSeg". The image(s) (polyT, DAPI, and cell boundaries)
#' are also read as \code{\link{SpatRaster}} objects so they are not loaded into
#' memory unless necessary. Because the image's origin is the top left while the
#' geometry's origin is bottom left, either the image or the geometry needs to
#' be flipped. Because the image accompanying MERFISH datasets are usually very
#' large, the coordinates will be flipped so the flipping operation won't load
#' the entire image into memory. Large datasets with hundreds of thousands of
#' cells can take a while to read if reading transcript spots as it takes a
#' while to convert the spots to MULTIPOINT geometries.
#'
#' @inheritParams SpatialFeatureExperiment
#' @inheritParams formatTxSpots
#' @param data_dir Top level directory of Vizgen output, which contains
#'   directories \code{cell_boundaries} and \code{images}, and files
#'   \code{cell_by_gene.csv}, \code{cell_metadata.csv}, and
#'   \code{detected_transcripts.csv}.
#' @param z Integer, z index to read, only affecting images read since cell
#'   segmentation for all z-planes are the same and cell centroids are only
#'   provided in 2 dimensions.
#' @param max_flip Maximum size of the image allowed to flip the image. Because
#'   the image will be loaded into memory to be flipped. If the image is larger
#'   than this size then the coordinates will be flipped instead.
#' @param flip To flip the image, geometry coordinates, or none. Because the
#'   image has the origin at the top left while the geometry has origin at the
#'   bottom left, one of them needs to be flipped for them to match. If one of
#'   them is already flipped, then use "none". The image will not be flipped if
#'   it's GeoTIFF.
#' @param image Which image(s) to load, can be "DAPI", "PolyT", "Cellbound" or
#'   any combination of them.
#' @param min_area Minimum cell area in square microns. Anything smaller will be
#'   considered artifact or debris and removed.
#' @param filter_counts Logical, whether to keep cells with counts \code{> 0}.
#' @param add_molecules Logical, whether to add transcripts coordinates to an
#'   object.
#' @param use_bboxes If no segmentation output is present, use
#'   \code{cell_metadata} to make bounding boxes instead.
#' @param use_cellpose Whether to read the parquet files from CellPose cell
#'   segmentation. If \code{FALSE}, cell segmentation will be read from the HDF5
#'   files. Note that reading HDF5 files for numerous FOVs is very slow.
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying parallel
#'   processing backend and number of threads to use for parallelizable tasks:
#'   \enumerate{ \item To load cell segmentation from HDF5 files from different
#'   fields of view (FOVs) with multiple cores. A progress bar can be configured
#'   in the \code{\link{BiocParallelParam}} object. When there are numerous
#'   FOVs, reading in the geometries can be time consuming, so we recommend
#'   using a server and larger number of threads. This argument is not used if
#'   \code{use_cellpose = TRUE} and the parquet file is present.
#'
#'   \item To get the largest piece and see if it's larger than \code{min_area}
#'   when there are multiple pieces in the cell segmentation for one cell.}
#' @param ... Other arguments passed to \code{\link{formatTxSpots}} to format
#'   and add the transcript spots data to the SFE object, except that extent is
#'   read from `manifest.json` and that `dest = "rowGeometry"` because the spot
#'   coordinates are in micron space and are not discrete so converting the
#'   transcript spots to raster won't work. A default is set for `file_out` to
#'   save the reformatted transcript spots to disk by default since reloading
#'   the reformatted form is much more efficient. Reading the original detected
#'   transcripts csv file can take up a lot of memory. Expect at least twice the
#'   size of that csv file, even more if using multiple threads. So we STRONGLY
#'   recommend saving the reformatted results to disk.
#' @concept Read data into SFE
#' @return A \code{SpatialFeatureExperiment} object.
#' @export
#' @note Since the transcript spots file is often very large, we recommend only
#'   using \code{add_molecules = TRUE} on servers with a lot of memory. If
#'   reading all z-planes, conversion of transcript spot geometry to parquet
#'   file might fail. In this case, call \code{formatTxSpots} separately, use
#'   \code{split_genes = TRUE} to save spots of each gene to a separate parquet
#'   file, and later selectively load genes of interest.
#' @importFrom sf st_area st_geometry<- st_as_sf st_read
#' @importFrom terra rast ext vect
#' @importFrom BiocParallel bpmapply bplapply
#' @importFrom rlang check_installed
#' @importFrom SpatialExperiment imgData<-
#' @importFrom data.table fread merge.data.table rbindlist is.data.table
#' @examples
#' dir_use <- system.file("extdata/vizgen_cellbound", package = "SpatialFeatureExperiment")
#' sfe <- readVizgen(dir_use, z = 3L, image = "PolyT",
#' flip = "geometry")
#'
#' ## Filtering of counts, and addition of molecule coordinates..
#' sfe <- readVizgen(dir_use, z = 3L, image = "PolyT", filter_counts = TRUE,
#' add_molecules = TRUE, flip = "geometry")
readVizgen <- function(data_dir,
                       z = 3L,
                       sample_id = "sample01", # How often do people read in multiple samples?
                       min_area = 15,
                       image = c("DAPI", "PolyT", "Cellbound"),
                       flip = c("geometry", "image", "none"),
                       max_flip = "50 MB",
                       filter_counts = FALSE,
                       add_molecules = FALSE,
                       use_bboxes = FALSE,
                       use_cellpose = TRUE,
                       BPPARAM = SerialParam(),
                       file_out = file.path(data_dir, "detected_transcripts.parquet"), ...) {
    data_dir <- normalizePath(data_dir, mustWork = TRUE)
    flip <- match.arg(flip)
    image <- match.arg(image, several.ok = TRUE)
    if ((any(z < 0) || any(z > 6)) && z != "all") {
        stop("z must be beween 0 and 6 (inclusive).")
    }

    # Read images----------
    # sanity on image names
    # .."Cellbound" image usually has a digit, eg "Cellbound3"
    image_regex <- image
    if (any("Cellbound" %in% image)) {
        image_regex[which(image %in% "Cellbound")] <-
            paste0(grep("Cell", image_regex, value = TRUE), "\\d") }

    if (z == "all") {
        img_pattern <- paste0("mosaic_(", paste(image_regex, collapse = "|"), ")_z-?\\d+\\.tiff?$")
    } else {
        num_pattern <- paste(z, collapse = "|")
        img_pattern <- paste0("mosaic_(", paste(image_regex, collapse = "|"), ")_z",
                              num_pattern, "\\.tiff?$")
    }
    img_fn <- list.files(file.path(data_dir, "images"), pattern = img_pattern,
                         ignore.case = TRUE, full.names = TRUE)
    if_exists <- vapply(image, function(img) any(grepl(img, img_fn, ignore.case = TRUE)),
                        FUN.VALUE = logical(1))
    if (!all(if_exists)) {
        warning("The image file(s) for ", "`", paste0(image[!if_exists], collapse = "|"), "`",
                " in this z-plane don't exist, or have non-standard file name(s).")
    }
    do_flip <- .if_flip_img(img_fn, max_flip)
    if (!length(img_fn)) flip <- "none"
    else if (!any(do_flip) && flip == "image") flip <- "geometry"

    # Read cell segmentation-------------
    # Use segmentation output from ".parquet" file
    # check if ".parquet" file is present
    parq <-
        # look in the current directory
        list.files(data_dir,
                   pattern = ".parquet$",
                   full.names = TRUE)
    if (length(parq) == 0 || any(grepl("detected_transcripts", parq))) {
        # look in the sub directory (if nothing is found)
        parq <- list.files(data_dir,
                           pattern = ".parquet$",
                           full.names = TRUE,
                           recursive = TRUE)
    }

    # set to use .parquet" file if present
    use.parquet <- any(length(parq)) & use_cellpose
    if (use.parquet) {
        # sanity check
        parq_sanity <-
            grepl("cell_boundaries|micron_space", parq)
        # make sure only single file is read
        if (length(parq) >= 1 && any(parq_sanity)) {
            # eg, if two files are present:
            # `cellpose_micron_space.parquet`
            # `cellpose_mosaic_space.parquet`
            # or any other `parquet` files
            # use µm units
            parq_clean <-
                grep("cell_boundaries|micron_space",
                     parq, value = TRUE)
            message(">>> ", length(parq), " `.parquet` files exist:",
                    paste0("\n", parq))
            parq <- parq_clean
            if (any(grepl("cell_boundaries.parquet", parq))) {
                # use default segmentaion file
                parq <- grep("cell_boundaries.parquet", parq, value = TRUE)
            } else if (any(grepl("hdf5s_micron", parq))) {
                # use previously processed/saved `hdf5` files
                parq <- grep("hdf5s_micron", parq, value = TRUE) }
            # final sanity
            if (length(parq) > 1) {
                stop("only 1 `.parquet` file can be read, check `data_dir` content")
            }
            message(">>> using -> " , parq)
        } else if (all(parq_sanity == FALSE)) { parq <- NULL }
        if (!is.null(parq)) {
            message(">>> Cell segmentations are found in `.parquet` file",
                    if (any(grepl("hdf5s_micron", parq))) {
                        paste0("\n", ">>> processed hdf5 files will be used") })
            fn <- parq
            # read file and filter to keep selected single z section as they're the same anyway
            polys <- st_read(fn, int64_as_string = TRUE, crs = NA, quiet = TRUE)
            # Can use any z-plane since they're all the same
            # This way so this part still works when the parquet file is written after
            # reading in HDF5 the first time. Only writing one z-plane to save disk space.
            polys <- polys[polys$ZIndex == polys$ZIndex[1],]
            polys$ZIndex <- if (z == "all") 3L else z[1]
            # filtering cell polygons
            polys <- .filter_polygons(polys, min_area,
                                      BPPARAM = BPPARAM)
            st_geometry(polys) <- "geometry"
            if ("EntityID" %in% names(polys))
                polys$ID <- polys$EntityID
            if (!"ZLevel" %in% names(polys)) # For reading what's written after HDF5
                polys$ZLevel <- 1.5 * (polys$ZIndex + 1L)
            polys <- polys[,c("ID", "ZIndex", "Type", "ZLevel", "geometry")]
        } else {
            warning("No '.parquet' or `hdf5` files present, check input directory -> `data_dir`")
            polys <- NULL }
    } else {
        check_installed("rhdf5")
        fns <- list.files(file.path(data_dir, "cell_boundaries"),
                          "*.hdf5", full.names = TRUE)
        if (length(fns)) {
            message(">>> Reading '.hdf5' files..")
            polys <- bpmapply(.h52poly_fov, fn = fns, SIMPLIFY = FALSE,
                              BPPARAM = BPPARAM,
                              # use mid z section
                              MoreArgs = list(z = ifelse(z == "all", 3L, z)))
            polys <- if (length(polys) == 1L) polys[[1]] else rbindlist(polys) |> st_as_sf()
            # recalculate bbox since rbindlist() takes the 1st polygon's bbox
            # see this -> https://github.com/Rdatatable/data.table/issues/4681
            polys <- polys[1:nrow(polys),]
            polys$Type <- "cell"
            parq_file <- file.path(data_dir, "hdf5s_micron_space.parquet")
            if (!file.exists(parq_file)) {
                st_write(polys, dsn = parq_file, driver = "Parquet")
            }
        } else if (length(fns) == 0) {
            warning("No '.hdf5' files present, check input directory -> `data_dir`")
            polys <- NULL
        }
    }
    if (!is.null(polys) && nrow(polys) == 0L)
        stop("No polygons left after filtering.")
    if (flip == "geometry" && !is.null(polys)) {
        # Flip the coordinates
        mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
        st_geometry(polys) <- st_geometry(polys) * mat_flip
    }

    # get count data file
    mat_fn <- .check_vizgen_fns(data_dir, "cell_by_gene")

    # Column without colname is read as V1
    mat <- fread(mat_fn, colClasses = list(character = 1))
    # TODO: write the sparse matrix to disk as hdf5 or mtx and check for the file
    # next time, so the costly step using a lot of memory to read in the dense
    # matrix can be avoided
    # get spatial metadata file---------
    meta_fn <- .check_vizgen_fns(data_dir, "cell_metadata")
    metadata <- fread(meta_fn, colClasses = list(character = 1))
    if (any(names(metadata) == "transcript_count") && filter_counts) {
        message(">>> ..filtering `cell_metadata` - keep cells with `transcript_count` > 0")
        metadata <- metadata[metadata$transcript_count > 0,]
    }

    if (!is.null(polys)) {
        # remove NAs when matching
        metadata <-
            metadata[match(polys$ID, metadata[[1]]) |> na.omit(),]
    }
    rownames(metadata) <- metadata[[1]]
    metadata[,1] <- NULL
    if (flip == "geometry") {
        metadata$center_y <- -metadata$center_y
    }

    # convert counts df to sparse matrix------------
    mat <- mat[match(rownames(metadata), mat[[1]]),] # polys already matched to metadata
    rns <- mat[[1]]
    mat[,1] <- NULL
    mat <- mat |>
        as.matrix() |>
        as("CsparseMatrix") |> # convert to sparse matrix
        Matrix::t() # transpose sparse matrix
    colnames(mat) <- rns
    # If metadata isn't already filtered
    if (!"transcript_count" %in% names(metadata) && filter_counts) {
        inds <- colSums(mat) > 0
        mat <- mat[,inds]
        metadata <- metadata[inds,]
        polys <- polys[inds,]
    }

    # check matching cell ids in polygon geometries, should match the count matrix cell ids
    if (!is.null(polys) &&
        !identical(polys$ID, rns)) {
        # filter geometries
        matched.cells <- match(rns, polys$ID) |> na.omit()
        message(">>> filtering geometries to match ", length(matched.cells),
                " cells with counts > 0")
        polys <- polys[matched.cells, , drop = FALSE]
    }

    if (any(if_exists)) {
        manifest <- fromJSON(file = file.path(data_dir, "images", "manifest.json"))
        extent <- setNames(manifest$bbox_microns, c("xmin", "ymin", "xmax", "ymax"))
        if (flip == "geometry") {
            extent[c("ymin", "ymax")] <- -extent[c("ymax", "ymin")]
        }
        # Set up ImgData
        img_dfs <- lapply(img_fn, function(fn) {
            # Now allowing multiple z planes
            id_use <- sub("^mosaic_", "", basename(fn))
            id_use <- sub("\\.tiff?$", "", id_use)
            .get_imgData(fn, sample_id = sample_id, image_id = id_use,
                         extent = extent, flip = (flip == "image"))
        })
        img_df <- do.call(rbind, img_dfs)
    }
    sfe <- SpatialFeatureExperiment(assays = list(counts = mat),
                                    colData = metadata,
                                    sample_id = sample_id,
                                    spatialCoordsNames = c("center_x", "center_y"),
                                    unit = "micron")

    # If none of segmentations are present, make bounding boxes
    # NOTE: might take some time to run
    if (use_bboxes && is.null(polys)) {
        message(">>> Creating bounding boxes from `cell_metadata`")
        # TODO: rewrite to use the now much faster df2sf ----
        bboxes_sfc <-
            bplapply(seq_len(nrow(metadata)),
                     function(i) {
                         bounds <- metadata[i, c("min_x", "max_x", "min_y", "max_y")]
                         names(bounds) <- c("xmin", "xmax", "ymin", "ymax")
                         st_as_sfc(st_bbox(bounds))
                     }, BPPARAM = BPPARAM)
        bboxes <- st_sf(geometry = st_sfc(bboxes_sfc))
        rownames(bboxes) <- rownames(metadata)
        cellSeg(sfe) <- bboxes
    }

    if (!is.null(polys)) {
        # sanity on geometries
        polys <- .check_st_valid(polys)
        rownames(polys) <- polys$ID
        polys$ID <- NULL
        cellSeg(sfe) <- polys
    }

    if (any(if_exists)) { imgData(sfe) <- img_df }

    if (add_molecules) {
        message(">>> Reading transcript coordinates")
        # get molecule coordiantes file
        mols_fn <- .check_vizgen_fns(data_dir, "detected_transcripts.csv")
        sfe <- addTxSpots(sfe, mols_fn, sample_id, BPPARAM = BPPARAM, z = z,
                          file_out = file_out, flip = (flip == "geometry"), ...)
    }
    sfe
}

.mols2geo <- function(mols, dest, spatialCoordsNames, gene_col, cell_col, not_in_cell_id) {
    # For one part of the split, e.g. cell compartment
    if (dest == "rowGeometry") {
        # Should have genes as row names
        # RAM concerns for parallel processing, wish I can stream
        mols <- df2sf(mols, geometryType = "MULTIPOINT",
                      spatialCoordsNames = spatialCoordsNames,
                      group_col = gene_col)
    } else {
        mols <- split(mols, mols[[gene_col]])
        mols <- lapply(mols, df2sf, geometryType = "MULTIPOINT",
                       spatialCoordsNames = spatialCoordsNames,
                       group_col = cell_col)
        names(mols) <- paste(names(mols), "spots", sep = "_")
    }
    mols
}

.mols2geo_split <- function(mols, dest, spatialCoordsNames, gene_col, cell_col,
                            not_in_cell_id, split_col) {
    if (!is.null(split_col) && split_col %in% names(mols)) {
        mols <- split(mols, mols[[split_col]])
        mols <- lapply(mols, .mols2geo, dest = dest,
                       spatialCoordsNames = spatialCoordsNames,
                       gene_col = gene_col, cell_col = cell_col,
                       not_in_cell_id = not_in_cell_id)
        if (dest == "colGeometry") {
            # Will be a nested list
            mols <- unlist(mols, recursive = FALSE)
            # names will be something like nucleus.Gapdh if split by compartment
        }
    } else {
        mols <- .mols2geo(mols, dest, spatialCoordsNames, gene_col, cell_col,
                          not_in_cell_id)
    }
    mols
}

# helper function to convert from raw bytes to character
.rawToChar_df <- function(input_df, BPPARAM = SerialParam()) {
    convert_ids <-
        lapply(input_df, function(x) is(x, "arrow_binary")) |> unlist() |> which()
    if (any(convert_ids)) {
        message(">>> Converting columns with raw bytes (ie 'arrow_binary') to character")
        cols_converted <-
            lapply(seq(convert_ids), function(i) {
                bplapply(input_df[,convert_ids][[i]], function(x) {
                    x <- rawToChar(x)
                }, BPPARAM = BPPARAM)
            })
        # replace the converted cell ids
        for (i in seq(cols_converted)) {
            input_df[,convert_ids][[i]] <- unlist(cols_converted[[i]])
        }
    }
    if (!is(input_df, "data.table")) {
        input_df <- data.table::as.data.table(input_df)
    }
    return(input_df)
}

.make_sql_query <- function(fn, gene_select, gene_col) {
    if (is.null(gene_select)) return(NA)
    lyr_name <- basename(fn) |> file_path_sans_ext()
    gene_part <- paste0("('", paste0(gene_select, collapse = "','"), "')")
    paste0("SELECT * FROM ", lyr_name, " WHERE ", gene_col, " IN ", gene_part)
}

.read_tx_output <- function(file_out, z, z_option, gene_col, return,
                            gene_select = NULL) {
    file_out <- normalizePath(file_out, mustWork = FALSE)
    file_dir <- file_path_sans_ext(file_out)
    # File or dir already exists, skip processing
    # read transcripts from ./detected_transcripts
    if (dir.exists(file_dir)) {
        # Multiple files
        pattern <- "\\.parquet$"
        # Need to deal with z-planes
        if (z != "all") {
            pattern <- paste0("_z", paste0(z, collapse = "|"), pattern)
        }
        fns <- list.files(file_dir, pattern, full.names = TRUE)
        if (!length(fns) && (length(z) == 1L || z_option == "3d")) {
            pattern <- "\\.parquet$"
            fns <- list.files(file_dir, pattern, full.names = TRUE)
        }
        if (length(fns)) {
            if (!return) return(file_dir)
            out <- lapply(fns, function(x) {
                q <- .make_sql_query(x, gene_select, gene_col)
                out <- st_read(x, query = q, int64_as_string = TRUE, quiet = TRUE,
                               crs = NA)
                out
            })
            # add names to a list
            names(out) <- gsub(".parquet", "",
                               x = list.files(file_dir, pattern))
            out <- lapply(out, function(x) {
                # row names are dropped in st_read/write_parquet
                rownames(x) <- x[[gene_col]]
                return(x)
            })
            return(out)
        }
        # read transcripts from detected_transcripts.parquet
    } else if (file.exists(file_out) && !dir.exists(file_dir)) {
        if (!return) return(file_out)
        out <- st_read(file_out, query = .make_sql_query(file_out, gene_select, gene_col),
                       int64_as_string = TRUE, quiet = TRUE, crs = NA)
        rownames(out) <- out[[gene_col]]
        return(out)
    }
}

#' Read and process transcript spots geometry for SFE
#'
#' The function `formatTxSpots` reads the transcript spot coordinates of
#' smFISH-based data and formats the data. The data is not added to an SFE
#' object. If the file specified in `file_out` already exists, then this file
#' will be read instead of the original file in the `file` argument, so the
#' processing is not run multiple times. The function `addTxSpots` adds the data
#' read and processed in `formatTxSpots` to the SFE object, and reads all
#' transcript spot data. To only read a subset of transcript spot data, first
#' use `formatTxSpots` to write the re-formatted data to disk. Then read the
#' specific subset and add them separately to the SFE object with the setter
#' functions.
#'
#' @param sfe A `SpatialFeatureExperiment` object.
#' @param file File with the transcript spot coordinates. Should be one row per
#'   spot when read into R and should have columns for coordinates on each axis,
#'   gene the transcript is assigned to, and optionally cell the transcript is
#'   assigned to. Must be csv, tsv, or parquet.
#' @param dest Where in the SFE object to store the spot geometries. This
#'   affects how the data is processed. Options: \describe{
#'   \item{rowGeometry}{All spots for each gene will be a `MULTIPOINT` geometry,
#'   regardless of whether they are in cells or which cells they are assigned
#'   to.} \item{colGeometry}{The spots for each gene assigned to a cell of
#'   interest will be a `MULTIPOINT` geometry; since the gene count matrix is
#'   sparse, the geometries are NOT returned to memory.}}
#' @param spatialCoordsNames Column names for the x, y, and optionally z
#'   coordinates of the spots. The defaults are for Vizgen.
#' @param z Index of z plane to read. Can be "all" to read all z-planes into
#'   MULTIPOINT geometries with XYZ coordinates. If z values are not integer,
#'   then spots with all z values will be read.
#' @param gene_col Column name for genes.
#' @param cell_col Column name for cell IDs, ignored if `dest = "rowGeometry"`.
#'   Can have length > 1 when multiple columns are needed to uniquely identify
#'   cells, in which case the contents of the columns will be concatenated, such
#'   as in CosMX data where cell ID is only unique within the same FOV. Default
#'   "cell_id" is for Vizgen MERFISH. Should be `c("cell_ID", "fov")` for CosMX.
#' @param not_in_cell_id Value of cell ID indicating that the spot is not
#'   assigned to any cell, such as "-1" in Vizgen MERFISH and "0" in CosMX. When
#'   there're multiple columns for `cell_col`, the first column is used to
#'   identify spots that are not in cells.
#' @param phred_col Column name for Phred scores of the spots.
#' @param min_phred Minimum Phred score to keep spot. By default 20, the
#'   conventional threshold indicating "acceptable", meaning that there's 1%
#'   chance that the spot was decoded in error.
#' @param split_col Categorical column to split the geometries, such as cell
#'   compartment the spots are assigned to as in the "CellComp" column in CosMX
#'   output.
#' @param file_out Name of file to save the geometry or raster to disk.
#'   Especially when the geometries are so large that it's unwieldy to load
#'   everything into memory. If this file (or directory for multiple files)
#'   already exists, then the existing file(s) will be read, skipping the
#'   processing. When writing the file, extensions supplied are ignored and
#'   extensions are determined based on `dest`.
#' @param return Logical, whether to return the geometries in memory. This does
#'   not depend on whether the geometries are written to file. Always `FALSE`
#'   when `dest = "colGeometry"`.
#' @param flip Logical, whether to flip the geometry to match image. Here the y
#'   coordinates are simply set to -y, so the original bounding box is not
#'   preserved. This is consistent with \code{readVizgen} and \code{readXenium}.
#' @param z_option What to do with z coordinates. "3d" is to construct 3D
#'   geometries. "split" is to create a separate 2D geometry for each z-plane so
#'   geometric operations are fully supported but some data wrangling is
#'   required to perform 3D analyses. When the z coordinates are not integers,
#'   3D geometries will always be constructed since there are no z-planes to
#'   speak of. This argument does not apply when `spatialCoordsNames` has length
#'   2.
#' @param BPPARAM \code{\link{BiocParallelParam}} object to specify
#'   multithreading to convert raw char in some parquet files to R objects. Not
#'   used otherwise.
#' @param sample_id Which sample in the SFE object the transcript spots should
#'   be added to.
#' @param gene_select Character vector of a subset of genes as appearing in the
#'   column \code{gene_col} and consistent with the row names of \code{sfe} to
#'   add. If \code{NULL}, then all genes that have transcript spots are added.
#'   Only relevant when reading data from formatted files on disk.
#' @return A sf data frame for vector geometries if `file_out` is not set.
#'   `SpatRaster` for raster. If there are multiple files written, such as when
#'   splitting by cell compartment or when `dest = "colGeometry"`, then a
#'   directory with the same name as `file_out` will be created (but without the
#'   extension) and the files are written to that directory with informative
#'   names. `parquet` files that can be read with `st_read` is written for
#'   vector geometries and `geotiff` files that can be read with `terra::rast`
#'   is written for raster. When `return = FALSE`, the file name or directory
#'   (when there're multiple files) is returned.
#' @note When `dest = "colGeometry"`, the geometries are always written to disk
#'   and not returned in memory, because this is essentially the gene count
#'   matrix, which is sparse. This kind of reformatting is implemented so users
#'   can read in MULTIPOINT geometries with transcript spots for each gene
#'   assigned to each cell for spatial point process analyses, where not all
#'   genes are loaded at once.
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom terra nlyr
#' @export
#' @return The `sf` data frame, or path to file where geometries are written if
#'   `return = FALSE`.
#' @rdname formatTxSpots
#' @examples
#' # Default arguments are for MERFISH
#' dir_use <- system.file("extdata/vizgen_cellbound", package = "SpatialFeatureExperiment")
#' g <- formatTxSpots(file.path(dir_use, "detected_transcripts.csv"))
#'
#' # For CosMX, note the colnames, also dest = "colGeometry"
#' # Results are written to the tx_spots directory
#' dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
#' cg <- formatTxSpots(file.path(dir_use, "Run5642_S3_Quarter_tx_file.csv"),
#' dest = "colGeometry", z = "all",
#' cell_col = c("cell_ID", "fov"),
#' gene_col = "target", not_in_cell_id = "0",
#' spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
#' file_out = "tx_spots")
formatTxSpots <- function(file, dest = c("rowGeometry", "colGeometry"),
                          spatialCoordsNames = c("global_x", "global_y", "global_z"),
                          gene_col = "gene", cell_col = "cell_id", z = "all",
                          phred_col = "qv", min_phred = 20, split_col = NULL,
                          not_in_cell_id = c("-1", "UNASSIGNED"),
                          z_option = c("3d", "split"), flip = FALSE,
                          file_out = NULL, gene_select = NULL, BPPARAM = SerialParam(),
                          return = TRUE) {
    file <- normalizePath(file, mustWork = TRUE)
    dest <- match.arg(dest)
    z_option <- match.arg(z_option)
    ext <- file_ext(file)
    if (dest == "colGeometry") {
        return <- FALSE
        if (is.null(file_out))
            stop("file_out must be specified for dest = 'colGeometry'.")
    }
    if (!ext %in% c("csv", "gz", "tsv", "txt", "parquet")) {
        stop("The file must be one of csv, gz, tsv, txt, or parquet")
    }
    if (!is.null(file_out) && (file.exists(file_out) ||
        dir.exists(basename(file_out) |> file_path_sans_ext()))) {
        out <- .read_tx_output(file_out, z, z_option, gene_col, return)
        return(out)
    }
    if (!is.numeric(z) && z != "all") {
        stop("z must either be numeric or be 'all' indicating all z-planes.")
    }
    if (ext == "parquet") {
        check_installed("arrow")
        mols <- arrow::read_parquet(file)
        # convert cols with raw bytes to character
        # NOTE: can take a while.
        mols <- .rawToChar_df(mols, BPPARAM = BPPARAM)
        # sanity, convert to data.table
        if (!is(mols, "data.table")) {
            mols <- data.table::as.data.table(mols)
        }
    } else {
        mols <- fread(file)
    }
    ind <- !spatialCoordsNames[1:2] %in% names(mols)
    if (any(ind)) {
        col_offending <- setdiff(spatialCoordsNames[1:2], names(mols))
        ax <- c("x", "y")
        stop(paste(ax[ind], collapse = ", "), " coordinate column(s) ",
             paste(col_offending, collapse = ", "), " not found.")
    }
    spatialCoordsNames <- intersect(spatialCoordsNames, names(mols))
    if (flip) {
        y_name <- spatialCoordsNames[2]
        if (!is.data.table(mols)) ..y_name <- y_name
        mols[,y_name] <- -mols[,..y_name]
    }
    # Check z
    use_z <- length(spatialCoordsNames) == 3L
    if (use_z) {
        zs <- mols[[spatialCoordsNames[3]]]
        if (is.null(zs)) { # z column not found
            spatialCoordsNames <- spatialCoordsNames[-3]
            use_z <- FALSE
        }
        if (all(floor(zs) == zs)) { # integer z values
            if (z != "all") {
                if (all(!z %in% unique(zs)))
                    stop("z plane(s) specified not found.")
                inds <- mols[[spatialCoordsNames[3]]] %in% z
                mols <- mols[inds,, drop = FALSE]
                if (length(z) == 1L) {
                    spatialCoordsNames <- spatialCoordsNames[-3]
                    use_z <- FALSE
                }
            }
        } else {
            z <- "all" # Non-integer z values
            z_option <- "3d"
        }
    }
    if (phred_col %in% names(mols) && !is.null(min_phred)) {
        mols <- mols[mols[[phred_col]] >= min_phred,]
    }
    message(">>> Converting transcript spots to geometry")
    if (dest == "colGeometry") {
        if (!length(cell_col) || any(!cell_col %in% names(mols)))
            stop("Column indicating cell ID not found.")
        mols <- mols[!mols[[cell_col[1]]] %in% not_in_cell_id,]
        if (length(cell_col) > 1L) {
            if (!is.data.table(mols)) ..cell_col <- cell_col
            cell_col_use <- do.call(paste, c(mols[,..cell_col], sep = "_"))
            mols$cell_id_ <- cell_col_use
            mols[,cell_col] <- NULL
            cell_col <- "cell_id_"
        }
        mols[,cell_col] <- as.character(mols[[cell_col]])
    }
    if (!is.null(file_out)) {
        file_out <- normalizePath(file_out, mustWork = FALSE)
        file_dir <- file_path_sans_ext(file_out)
    }
    if (z_option == "split" && use_z) {
        mols <- split(mols, mols[[spatialCoordsNames[3]]])
        mols <- lapply(mols, .mols2geo_split, dest = dest,
                       spatialCoordsNames = spatialCoordsNames[1:2],
                       gene_col = gene_col, cell_col = cell_col,
                       not_in_cell_id = not_in_cell_id, split_col = split_col)
        # If list of list, i.e. colGeometry, or do split
        if (!is(mols[[1]], "sf")) {
            names_use <- lapply(names(mols), function(n) {
                names_int <- names(mols[[n]])
                paste0(names_int, "_z", n)
            }) |> unlist()
            mols <- unlist(mols, recursive = FALSE)
            names(mols) <- names_use
        } else if (!is.null(file_out)) {
            names(mols) <- paste0(basename(file_dir), "_z", names(mols))
        } else {
            names(mols) <-
                file_path_sans_ext(file) |>
                basename() |>
                paste0("_z", names(mols))
        }
    } else {
        mols <- .mols2geo_split(mols, dest, spatialCoordsNames, gene_col, cell_col,
                                not_in_cell_id, split_col)
    }

    if (!is.null(file_out)) {
        message(">>> Writing reformatted transcript spots to disk")
        if (!dir.exists(dirname(file_out)))
            dir.create(dirname(file_out))
        if (is(mols, "sf")) {
            st_write(mols, file_out, driver = "Parquet", quiet = TRUE)
            if (!return) return(file_out)
        } else {
            if (!dir.exists(file_dir)) dir.create(file_dir)
            suppressWarnings({
                bplapply(names(mols), function(n) {
                    name_use <- gsub("/", ".", n)
                    st_write(mols[[n]], file.path(file_dir, paste0(name_use, ".parquet")),
                             driver = "Parquet")
                }, BPPARAM = SerialParam(progressbar = TRUE))
            })
            if (!return) return(file_dir)
        }
    }
    return(mols)
}

#' @rdname formatTxSpots
#' @export
addTxSpots <- function(sfe, file, sample_id = 1L, gene_select = NULL,
                       spatialCoordsNames = c("global_x", "global_y", "global_z"),
                       gene_col = "gene", z = "all",
                       phred_col = "qv", min_phred = 20, split_col = NULL,
                       z_option = c("3d", "split"), flip = FALSE,
                       file_out = NULL, BPPARAM = SerialParam()) {
    sample_id <- .check_sample_id(sfe, sample_id)
    z_option <- match.arg(z_option)
    dest <- "rowGeometry"
    gene_select <- intersect(gene_select, rownames(sfe))
    if (!length(gene_select)) gene_select <- NULL
    if (!is.null(gene_select)) {
        if (is.null(file_out))
            stop("file_out must be specified if only loading a subset of genes")
        return <- FALSE
        partial <- TRUE
    } else {
        return <- TRUE
        partial <- FALSE
    }
    mols <- formatTxSpots(file, dest = dest, spatialCoordsNames = spatialCoordsNames,
                          gene_col = gene_col, z = z, phred_col = phred_col,
                          min_phred = min_phred, split_col = split_col,
                          flip = flip, gene_select = gene_select,
                          z_option = z_option, file_out = file_out,
                          BPPARAM = BPPARAM, return = return)
    if (is(mols, "sf")) {
        txSpots(sfe, withDimnames = TRUE, partial = partial) <- mols
    } else if (is(mols, "list")) {
        rowGeometries(sfe, partial = partial) <- mols
    }

    # make sure that sfe and rowGeometries have the same features
    # NOTE, if `min_phred = NULL`, no filtering of features occur
    if (!is.null(min_phred)) {
        if (length(rowGeometries(sfe)) > 1) {
            # check if all features match between rowGeometries and SFE object
            gene_names <-
                lapply(rowGeometries(sfe), function(i) {
                    gene_indx <-
                        which(rownames(sfe) %in% na.omit(i[[gene_col]]))
                    gene_name <- rownames(sfe[gene_indx,])
                    return(gene_name)
                }) |> unlist() |> unique()
            if (!all(rownames(sfe) %in% gene_names)) {
                gene_indx <-
                    which(rownames(sfe) %in% gene_names)
                # genes/features that are removed
                genes_rm <- rownames(sfe)[-gene_indx]
                message(">>> Total of ", length(genes_rm),
                        " features/genes with no transcript detected or `min_phred` < ",
                        min_phred, " are removed from SFE object",
                        "\n", ">>> To keep all features -> set `min_phred = NULL`")
                sfe <- sfe[gene_indx,]
            }
        } else if (length(rowGeometries(sfe)) == 1) {
            # NOTE, transcripts are filtered with default qv/min_phred >= 20
            # for txSpots
            if (!all(rownames(sfe) %in% rowGeometry(sfe)[[gene_col]])) {
                # match gene names from rowGeometry removing NAs
                gene_indx <-
                    which(rownames(sfe) %in% rowGeometry(sfe)[[gene_col]] |> na.omit())
                genes_rm <- rownames(sfe)[-gene_indx]
                message(">>> Total of ", length(genes_rm),
                        " features/genes with `min_phred` < ", min_phred, " are removed from SFE object",
                        "\n", ">>> To keep all features -> set `min_phred = NULL`")
                # subset sfe to keep genes present in rowGeometry
                sfe <- sfe[gene_indx,]
            }
        }
    }
    sfe
}


#' Read CosMX data into SFE
#'
#' This function reads the standard CosMX output into an SFE object, as in
#' "Basic Data Files" on the Nanostring website.
#'
#' @inheritParams readVizgen
#' @param split_cell_comps Logical, whether to split transcript spot geometries
#'   by cell compartment. Only relevant when `add_molecules = TRUE`.
#' @return An SFE object. Cell polygons are written to
#'   `cell_boundaries_sf.parquet` in `data_dir`. If reading transcript spots
#'   (`add_molecules = TRUE`), then the reformatted transcript spots are saved
#'   to file specified in the `file_out` argument, which is by default
#'   `tx_spots.parquet` in the same directory as the rest of the data.
#' @export
#' @examples
#' dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
#' file.copy(dir_use, ".", recursive = TRUE)
#' sfe <- readCosMX("cosmx", z = "all", add_molecules = TRUE)
#' # Clean up
#' unlink("cosmx", recursive = TRUE)
readCosMX <- function(data_dir,
                      z = 3L,
                      sample_id = "sample01", # How often do people read in multiple samples?
                      add_molecules = FALSE,
                      split_cell_comps = FALSE,
                      BPPARAM = SerialParam(),
                      file_out = file.path(data_dir, "tx_spots.parquet"), ...) {
    data_dir <- normalizePath(data_dir, mustWork = TRUE)
    fns <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
    fn_metadata <- grep("metadata", fns, value = TRUE)
    fn_mat <- grep("exprMat", fns, value = TRUE)
    fn_polys <- grep("polygons", fns, value = TRUE)

    meta <- fread(fn_metadata)
    mat <- fread(fn_mat) # TODO: write to h5 or mtx. Consult alabaster.sce
    polys <- fread(fn_polys)

    meta$cell_ID <- paste(meta$cell_ID, meta$fov, sep = "_")
    mat$cell_ID <- paste(mat$cell_ID, mat$fov, sep = "_")
    polys$cellID <- paste(polys$cellID, polys$fov, sep = "_")

    mat <- mat[match(meta$cell_ID, mat$cell_ID),]
    cell_ids <- mat$cell_ID
    mat <- mat[,3:ncol(mat)] |>
        as.matrix() |>
        as("CsparseMatrix") |> Matrix::t()
    colnames(mat) <- cell_ids

    poly_sf_fn <- file.path(data_dir, "cell_boundaries_sf.parquet")
    if (file.exists(poly_sf_fn)) {
        message(">>> File cell_boundaries_sf.parquet found")
        polys <- st_read(poly_sf_fn, int64_as_string = TRUE, quiet = TRUE, crs = NA)
        rownames(polys) <- polys$cellID
    } else {
        message(">>> Constructing cell polygons")
        polys <- df2sf(polys, spatialCoordsNames = c("x_global_px", "y_global_px"),
                       geometryType = "POLYGON",
                       id_col = "cellID")
        polys <- polys[match(meta$cell_ID, polys$cellID),]
        st_write(polys, poly_sf_fn, driver = "Parquet")
    }

    sfe <- SpatialFeatureExperiment(list(counts = mat), colData = meta,
                                    spatialCoordsNames = c("CenterX_global_px", "CenterY_global_px"),
                                    unit = "full_res_image_pixel")
    # sanity on geometries
    polys <- .check_st_valid(polys)
    cellSeg(sfe) <- polys

    if (add_molecules) {
        message(">>> Reading transcript coordinates")
        fn <- grep("tx_file.csv", fns, value = TRUE)
        split_col <- if (split_cell_comps) "CellComp" else NULL
        sfe <- addTxSpots(sfe, fn, spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
                          gene_col = "target", split_col = split_col,
                          file_out = file_out, z = z,
                          BPPARAM = BPPARAM, ...)
    }
    sfe
}

.check_xenium_fns <- function(data_dir, keyword, no_raw_bytes = FALSE) {
    fn_all <-
        list.files(data_dir,
                   pattern = keyword,
                   full.names = TRUE)
    if (any(grep(keyword, fn_all))) {
        # Priorities: 1. _sf.parquet, 2. csv, 3. parquet
        #..since .parquet has cols with raw bytes format in v1.3 or less
        fn <- grep("_sf.parquet", fn_all, value = TRUE)
        if (no_raw_bytes) {
            if (!length(fn)) fn <- grep(".parquet", fn_all, value = TRUE)
            if (!length(fn)) fn <- grep(".csv", fn_all, value = TRUE)
        } else {
            if (!length(fn)) fn <- grep(".csv", fn_all, value = TRUE)
            if (!length(fn)) fn <- grep(".parquet", fn_all, value = TRUE)
        }
    }
    if (!length(fn)) {
        stop("No `", keyword, "` file is available")
    }
    fn
}

#' Read 10X Xenium output as SpatialFeatureExperiment
#'
#' This function reads the standard 10X Xenium output into an SFE object.
#'
#' @inheritParams readVizgen
#' @param image Which image(s) to load, can be "morphology_mip",
#'   "morphology_focus" or both. Note that in Xenium Onboarding Analysis (XOA)
#'   v2, there is no longer "morphology_mip" and "morphology_focus" is a
#'   directory with 4 images corresponding to 4 channels: DAPI, "Cadherin", 18S,
#'   and Vimentin. So this argument is ignored for XOA v2.
#' @param segmentations Which segmentation outputs to read, can be "cell",
#'   "nucleus", or both.
#' @param row.names String specifying whether to use Ensembl IDs ("id") or gene
#'   symbols ("symbol") as row names. If using symbols, the Ensembl ID will be
#'   appended to disambiguate in case the same symbol corresponds to multiple
#'   Ensembl IDs. Always "symbol" if `add_molecules = TRUE` because only gene
#'   symbols are used in the transcript spot files.
#' @return An SFE object. If reading segmentations, the cell or nuclei
#'   segmentation will be saved to `cell_boundaries_sf.parquet` and
#'   `nucleus_boundaries_sf.parquet` respectively in `data.dir` so next time the
#'   boundaries can be read much more quickly. If reading transcript spots
#'   (`add_molecules = TRUE`), then the reformatted transcript spots are saved
#'   to file specified in the `file_out` argument, which is by default
#'   `tx_spots.parquet` in the same directory as the rest of the data. If images
#'   are present, then the images will be of the \code{BioFormatsImage} class
#'   and not loaded into memory until necessary in later operations.
#' @note Sometimes when reading images, you will see this error the first time:
#' 'java.lang.NullPointerException: Cannot invoke
#' "loci.formats.DimensionSwapper.setMetadataFiltered(boolean)" because
#' "RBioFormats.reader" is null'. Rerun the code and it should work the second
#' time.
#' @export
#'
#' @importFrom sf st_area st_geometry<- st_as_sf st_write
#' @importFrom terra rast ext vect
#' @importFrom BiocParallel bpmapply bplapply
#' @importFrom rlang check_installed
#' @importFrom SpatialExperiment imgData<-
#' @importFrom SingleCellExperiment counts
#' @importFrom data.table fread merge.data.table rbindlist is.data.table
#' @importFrom DropletUtils read10xCounts
#' @examples
#' # TODO: Example code for Xenium toy data
#'
#' # custom example run:
#' if (FALSE)
#' sfe <-
#'  readXenium(data_dir = data_dir,
#'  sample_id = "test_xenium",
#'  image = c("morphology_focus", "morphology_mip"),
#'  segmentations = c("cell", "nucleus"),
#'  flip = "geometry",
#'  filter_counts = TRUE,
#'  add_molecules = TRUE,
#'  file_out = NULL)
#'
readXenium <- function(data_dir,
                       sample_id = "sample01",
                       image = c("morphology_focus", "morphology_mip"),
                       segmentations = c("cell", "nucleus"),
                       row.names = c("id", "symbol"),
                       flip = c("geometry", "image", "none"),
                       max_flip = "50 MB",
                       filter_counts = FALSE,
                       add_molecules = FALSE,
                       BPPARAM = SerialParam(),
                       file_out = file.path(data_dir, "tx_spots.parquet"), ...) {
    data_dir <- normalizePath(data_dir, mustWork = TRUE)
    flip <- match.arg(flip)
    image <- match.arg(image, several.ok = TRUE)
    row.names <- match.arg(row.names)
    if (add_molecules) {
        message(">>> Must use gene symbols as row names when adding transcript spots.")
        row.names <- "symbol"
    }
    experiment <- fromJSON(file = file.path(data_dir, "experiment.xenium"),
                           simplify = TRUE)
    xoa_version <- experiment$analysis_sw_version
    major_version <- substr(xoa_version, 8, 8) |> as.integer()
    minor_version <- substr(xoa_version, 10, 10) |> as.integer()

    # Read images-----------
    # supports 2 images, in XOA v1:
    # `morphology_mip.ome.tif` - 2D maximum projection intensity (MIP) image of the tissue morphology image.
    # `morphology_focus.ome.tif` - 2D autofocus projection image of the tissue morphology image.
    # XOA v2: morphology_focus directory with multi-file OME-TIFF
    if (major_version == 1L) {
        img_fn <-
            list.files(data_dir, full.names = TRUE,
                       pattern = "morphology_.*\\.ome\\.tif")
        if_exists <- vapply(image, function(img) any(grepl(img, img_fn, ignore.case = TRUE)),
                            FUN.VALUE = logical(1))
        if (!all(if_exists)) {
            warning("The image file(s) for ", "`", paste0(image[!if_exists], collapse = "|"), "`",
                    " don't exist, or have non-standard file name(s).")
        }
    } else { # For now there's only v2. We'll see what v3 will be like
        img_fn <- paste0("morphology_focus_000", 0:3, ".ome.tif")
        img_fn <- file.path(data_dir, "morphology_focus", img_fn)
        if_exists <- vapply(img_fn, file.exists, FUN.VALUE = logical(1))
        if (!all(if_exists)) {
            warning("Image file(s) ", paste0(img_fn[!if_exists], collapse = ", "),
                    " not found")
        }
    }
    use_imgs <- any(if_exists)
    do_flip <- .if_flip_img(img_fn, max_flip)
    if (!length(img_fn)) {
        flip <- "none"
    } else if (!any(do_flip) && flip == "image") { flip <- "geometry" }
    if (use_imgs) {
        # Set up ImgData
        if (major_version == 1L) {
            img_dfs <- lapply(img_fn, function(fn) {
                id_use <- sub("\\.ome\\.tiff?$", "", basename(fn))
                .get_imgData(fn, sample_id = sample_id,
                             image_id = id_use,
                             flip = (flip == "image"))
            })
            img_df <- do.call(rbind, img_dfs)
        } else {
            img_df <- .get_imgData(img_fn[1], sample_id = sample_id,
                                   image_id = "morphology_focus",
                                   flip = (flip == "image"))
        }
        if (flip == "geometry") {
            img_df$data <- lapply(img_df$data, function(img) {
                extent <- ext(img)
                img <- translateImg(img, v = c(0, extent["ymin"] - extent["ymax"]))
                img
            })
        }
    }

    # Read cell/nucleus segmentation ----
    if (!is.null(segmentations)) {
        # get files .parquet or .csv
        # What if only cell or only nucleus is available
        no_raw_bytes <- (major_version == 1L && minor_version > 4L) || major_version == 2L
        fn_segs <- c(cell = .check_xenium_fns(data_dir, "cell_boundaries", no_raw_bytes),
                     nucleus = .check_xenium_fns(data_dir, "nucleus_boundaries", no_raw_bytes))
        fn_segs <- fn_segs[segmentations]
        if (length(fn_segs) == 0) {
            warning("No segmentation files are found, check input directory -> `data_dir`")
            polys <- NULL
        }
        if (any(grep("_sf.parquet", fn_segs))) {
            message(">>> Preprocessed sf segmentations found\n",
                    ">>> Reading ", paste0(names(fn_segs), collapse = " and "),
                    " segmentations")
            # add cell id to rownames
            polys <- lapply(fn_segs, st_read, quiet = TRUE, int64_as_string = TRUE,
                            crs = NA)
        } else {
            if (no_raw_bytes) {
                if (any(grep("..parquet", fn_segs))) {
                    check_installed("arrow")
                    message(">>> Cell segmentations are found in `.parquet` file(s)", "\n",
                            ">>> Reading ", paste0(names(fn_segs), collapse = " and "),
                            " segmentations")
                    polys <- lapply(fn_segs, arrow::read_parquet)
                } else if (any(grep(".csv", fn_segs))) {
                    message(">>> Cell segmentations are found in `.csv` file(s)", "\n",
                            ">>> Reading ", paste0(names(fn_segs), collapse = " and "),
                            " segmentations")
                    # read .csv data
                    polys <- lapply(fn_segs, fread)
                }
            } else {
                if (any(grep(".csv", fn_segs))) {
                    message(">>> Cell segmentations are found in `.csv` file(s)", "\n",
                            ">>> Reading ", paste0(names(fn_segs), collapse = " and "),
                            " segmentations")
                    # read .csv data
                    polys <- lapply(fn_segs, fread)
                } else if (any(grep("..parquet", fn_segs))) {
                    check_installed("arrow")
                    message(">>> Cell segmentations are found in `.parquet` file(s)", "\n",
                            ">>> Reading ", paste0(names(fn_segs), collapse = " and "),
                            " segmentations")
                    polys <- lapply(fn_segs, arrow::read_parquet)
                    # convert cell ids, from raw bytes to character
                    polys <- lapply(polys, function(x)
                        .rawToChar_df(x, BPPARAM = BPPARAM))
                }
            }
            # generate sf dataframe with geometries
            # Flip the coordinates
            if (flip == "geometry" && !is.null(polys)) {
                polys <- lapply(polys, function(p) {
                    p$vertex_y <- -p$vertex_y
                    p
                })
            }
            instrument_version <- experiment$instrument_sw_version
            if (major_version == 2L && instrument_version != "Development") {
                if ("nucleus" %in% names(polys)) {
                    message(">>> Making MULTIPOLYGON nuclei geometries")
                    polys[["nucleus"]] <- df2sf(polys[["nucleus"]],
                                                c("vertex_x", "vertex_y"),
                                                id_col = "label_id",
                                                group_col = "cell_id",
                                                geometryType = "MULTIPOLYGON")
                }
                if ("cell" %in% names(polys)) {
                    message(">>> Making POLYGON cell geometries")
                    polys[["cell"]] <- df2sf(polys[["cell"]],
                                             c("vertex_x", "vertex_y"),
                                             id_col = "cell_id",
                                             geometryType = "POLYGON")
                }
            } else {
                message(">>> Making POLYGON geometries")
                polys <-
                    lapply(polys, function(x) {
                        df2sf(x, c("vertex_x", "vertex_y"), id_col = "cell_id",
                              geometryType = "POLYGON") })
            }

            fn_out <- c(cell = "cell_boundaries_sf.parquet",
                        nucleus = "nucleus_boundaries_sf.parquet")
            fn_out <- fn_out[names(fn_segs)]
            fn_out <- file.path(data_dir, fn_out)
            message(">>> Saving geometries to parquet files")
            for (i in seq_along(polys)) {
                st_write(polys[[i]], fn_out[[i]], driver = "Parquet", quiet = TRUE)
            }
        }
        # add names to polys list
        names(polys) <- c(cell = "cellSeg", nucleus = "nucSeg")[names(fn_segs)]
    } else { polys <- NULL }

    # Read metadata ----
    fn_meta <- .check_xenium_fns(data_dir, "cells.")
    if (length(fn_meta) == 0) {
        warning("No metadata files are found, check input directory -> `data_dir`")
        metadata <- NULL
    }
    if (no_raw_bytes) {
        if (any(grep(".parquet", fn_meta))) {
            check_installed("arrow")
            metadata <- arrow::read_parquet(fn_meta)
            message(">>> Reading cell metadata -> `cells.parquet`")
        } else if (any(grep(".csv", fn_meta))) {
            message(">>> Reading cell metadata -> `cells.csv`")
            # read .csv data
            metadata <- fread(fn_meta)
        }
    } else {
        if (any(grep(".csv", fn_meta))) {
            message(">>> Reading cell metadata -> `cells.csv`")
            # read .csv data
            metadata <- fread(fn_meta)
        } else if (any(grep(".parquet", fn_meta))) {
            check_installed("arrow")
            metadata <- arrow::read_parquet(fn_meta)
            message(">>> Reading cell metadata -> `cells.parquet`")
            # convert cell ids, from raw bytes to character
            metadata <- .rawToChar_df(metadata, BPPARAM = BPPARAM)
        }
    }

    # Read count matrix or SCE ----
    # all feature types are read in single count matrix and stored in rowData(sce)$Type
    #..ie -> 'Negative Control Probe, 'Negative Control Codeword', 'Unassigned Codeword'
    if (file.exists(file.path(data_dir, "cell_feature_matrix.h5"))) {
        message(">>> Reading h5 gene count matrix")
        sce <- read10xCounts(file.path(data_dir, "cell_feature_matrix.h5"),
                             col.names = TRUE, row.names = row.names)
    } else if (dir.exists(file.path(data_dir, "cell_feature_matrix"))) {
        message(">>> Reading mtx gene count matrix")
        sce <- read10xCounts(file.path(data_dir, "cell_feature_matrix"),
                             col.names = TRUE, row.names = row.names)
    } else { stop("No `cell_feature_matrix` files are found, check input directory -> `data_dir`") }

    # Filtering count matrix, metadata and segmentations ----
    # filtering metadata and count matrix
    if (any(names(metadata) == "transcript_counts") && filter_counts) {
        message(">>> ..filtering cell metadata - keep cells with `transcript_counts` > 0")
        metadata <- metadata[metadata$transcript_count > 0,]
        sce <- sce[,match(metadata$cell_id, colnames(sce)) |> na.omit()]
    } else {
        # if metadata isn't already filtered
        if (!"transcript_counts" %in% names(metadata) && filter_counts) {
            inds <- colSums(counts(sce)) > 0
            sce <- sce[,inds]
            metadata <- metadata[inds,]
        }}
    # filtering segmentations
    if (!is.null(polys)) {
        # polys should always be a list, even if it's length 1
        for (i in seq(polys)) {
            # filter geometries
            matched.cells <- match(colnames(sce), polys[[i]]$cell_id) |> na.omit()
            message(">>> filtering ", names(polys)[i],
                    " geometries to match ",
                    length(matched.cells), " cells with counts > 0")
            polys[[i]] <- polys[[i]][matched.cells, , drop = FALSE] }
    }
    metadata <- as.data.frame(metadata) |> as("DataFrame")
    rownames(metadata) <- metadata$cell_id
    metadata[,1] <- NULL
    if (flip == "geometry") {
        metadata$y_centroid <- -metadata$y_centroid
    }

    # Make SFE object ----
    colData(sce) <- metadata
    sfe <- toSpatialFeatureExperiment(sce, sample_id = sample_id,
                                      spatialCoordsNames = c("x_centroid", "y_centroid"),
                                      unit = "micron")

    # add segmentation geometries
    if (!is.null(polys)) {
        message(">>> Checking polygon validity")
        # sanity on geometries
        polys <- lapply(polys, .check_st_valid)
        polys <-
            lapply(polys, function(i) {
                rownames(i) <- i$cell_id
                i$cell_id <- NULL
                return(i)}
            )
        colGeometries(sfe) <- c(colGeometries(sfe), polys)
    }

    # add images
    if (use_imgs) imgData(sfe) <- img_df

    # Read transcript coordinates ----
    # NOTE z-planes are non-integer, cannot select or use `z` as in `readVizgen`
    if (add_molecules) {
        message(">>> Reading transcript coordinates")
        # get molecule coordiantes file
        fn_mols <- .check_xenium_fns(data_dir, "transcripts", no_raw_bytes)
        sfe <- addTxSpots(sfe, fn_mols,
                          sample_id,
                          gene_col = "feature_name",
                          spatialCoordsNames = c("x_location", "y_location", "z_location"),
                          z_option = "3d",
                          BPPARAM = BPPARAM, flip = (flip == "geometry"),
                          file_out = file_out, ...)
    }
    sfe
}
