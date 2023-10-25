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
  if (length(empty.inds)) { warning(">>> ..removing ",
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
#' @param data_dir Top level directory of Vizgen output, which contains
#'   directories \code{cell_boundaries} and \code{images}, and files
#'   \code{cell_by_gene.csv}, \code{cell_metadata.csv}, and
#'   \code{detected_transcripts.csv}.
#' @param z Index of z plane to read. Can be "all" to read all z-planes into
#'   MULTIPOINT geometries with XYZ coordinates. If z values are not integer,
#'   then spots with all z values will be read.
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
#' @param filter_counts Keep cells with counts \code{> 0}.
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
#'   \item To filter cell segmentation polygons based on size.
#'
#'   \item To convert centroid points to sf POINT geometry.
#'
#'   \item When the cell segmentation is unavailable, to convert cell bounding
#'   boxes to POLYGON geometries.
#'
#'   \item Convert transcript coordinates into MULTIPOLYGON for rowGeometry. }
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
#'   using \code{add_molecules = TRUE} on servers with a lot of memory.
#' @importFrom sf st_area st_geometry<- st_as_sf
#' @importFrom terra rast ext vect
#' @importFrom BiocParallel bpmapply bplapply
#' @importFrom rlang check_installed
#' @importFrom SpatialExperiment imgData<-
#' @importFrom data.table fread merge.data.table rbindlist is.data.table
#' @examples
#' dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
#' sfe <- readVizgen(dir_use, z = 0L, image = "PolyT",
#' flip = "geometry")
#'
#' ## Filtering of counts, and addition of molecule coordinates..
#' sfe <- readVizgen(dir_use, z = 0L, image = "PolyT", filter_counts = TRUE,
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
    img_pattern <- paste0("mosaic_(", paste(image_regex, collapse = "|"), ")_z-?\\d+\\.tif$")
  } else {
    num_pattern <- paste(z, collapse = "|")
    img_pattern <- paste0("mosaic_(", paste(image_regex, collapse = "|"), ")_z",
                          num_pattern, "\\.tif$")
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
      message(">>> ", length(parq), " `.parquet` files exists:",
              paste0("\n", parq), "\n", ">>> using -> " , parq_clean)
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
    } else if (all(parq_sanity == FALSE)) { parq <- NULL }
    if (!is.null(parq)) {
      rlang::check_installed("sfarrow")
      message(">>> Cell segmentations are found in `.parquet` file",
              if (any(grepl("hdf5s_micron", parq))) {
                paste0("\n", ">>> processed hdf5 files will be used") })
      fn <- parq
      # read file and filter to keep selected single z section as they're the same anyway
      polys <- sfarrow::st_read_parquet(fn)
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
        warning("No '.parquet' files present, check input directory -> `data_dir`")
        polys <- NULL }
    } else {
    rlang::check_installed("rhdf5")
    fns <- list.files(file.path(data_dir, "cell_boundaries"),
                      "*.hdf5", full.names = TRUE)
    if (length(fns)) {
      message(">>> Reading '.hdf5' files..")
      polys <- bpmapply(.h52poly_fov, fn = fns, SIMPLIFY = FALSE,
                        BPPARAM = BPPARAM,
                        # use mid z section
                        MoreArgs = list(z = ifelse(z == "all", 3L, z)))
      polys <- if (length(polys) == 1L) polys[[1]] else rbindlist(polys) |> st_as_sf()
      polys$Type <- "cell"
      parq_file <- file.path(data_dir, "hdf5s_micron_space.parquet")
      if (!file.exists(parq_file)) {
        suppressWarnings(sfarrow::st_write_parquet(polys, dsn = parq_file))
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
      metadata[match(polys$ID, metadata[[1]]) |> stats::na.omit(),]
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
    matched.cells <- match(rns, polys$ID) |> stats::na.omit()
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
      id_use <- sub("\\.tif$", "", id_use)
      .get_imgData(fn, sample_id = sample_id, image_id = id_use,
                   extent = extent, flip = (flip == "image"))
    })
    img_df <- do.call(rbind, img_dfs)
  }
  sfe <- SpatialFeatureExperiment(assays = list(counts = mat),
                                  colData = metadata,
                                  sample_id = sample_id,
                                  spatialCoordsNames = c("center_x", "center_y"),
                                  unit = "micron", BPPARAM = BPPARAM)

  # If none of segmentations are present, make bounding boxes
  # NOTE: might take some time to run
  if (use_bboxes && is.null(polys)) {
    message(">>> Creating bounding boxes from `cell_metadata`")
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
                      file_out = file_out, ...)
  }
  sfe
}

.mols2geo <- function(mols, dest, spatialCoordsNames, gene_col, cell_col, BPPARAM, not_in_cell_id) {
  # For one part of the split, e.g. cell compartment
  if (dest == "rowGeometry") {
    # Should have genes as row names
    # RAM concerns for parallel processing, wish I can stream
    mols <- df2sf(mols, geometryType = "MULTIPOINT", BPPARAM = BPPARAM,
                  spatialCoordsNames = spatialCoordsNames,
                  group_col = gene_col)
  } else {
    mols <- split(mols, mols[[gene_col]])
    mols <- bplapply(mols, df2sf, geometryType = "MULTIPOINT",
                     spatialCoordsNames = spatialCoordsNames,
                     group_col = cell_col,
                     # Does it get passed to df2sf? I don't think so, so parallelizing over cells
                     BPPARAM = BPPARAM)
    names(mols) <- paste(names(mols), "spots", sep = "_")
  }
  mols
}

.mols2geo_split <- function(mols, dest, spatialCoordsNames, gene_col, cell_col,
                            BPPARAM, not_in_cell_id, split_col) {
    if (!is.null(split_col) && split_col %in% names(mols)) {
        mols <- split(mols, mols[[split_col]])
        mols <- lapply(mols, .mols2geo, dest = dest, BPPARAM = BPPARAM,
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
                          BPPARAM, not_in_cell_id)
    }
    mols
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
#' @inheritParams readVizgen
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
#' @param z_option What to do with z coordinates. "3d" is to construct 3D
#'   geometries. However, note that GEOS underpinning \code{sf} does not support
#'   3D geometric operations. "split" is to create a separate 2D geometry for
#'   each z-plane so geometric operations are fully supported but some data
#'   wrangling is required to perform 3D analyses, and this is preferred when z
#'   coordinates have arbitrary units different from those of x and y. When the
#'   z coordinates are not integers, 3D geometries will always be constructed
#'   since there are no z-planes to speak of. This argument does not apply when
#'   `spatialCoordsNames` has length 2.
#' @return A sf data frame for vector geometries if `file_out` is not set.
#'   `SpatRaster` for raster. If there are multiple files written, such as when
#'   splitting by cell compartment or when `dest = "colGeometry"`, then a
#'   directory with the same name as `file_out` will be created (but without the
#'   extension) and the files are written to that directory with informative
#'   names. `parquet` files that can be read with `sfarrow::st_read_parquet` is
#'   written for vector geometries and `geotiff` files that can be read with
#'   `terra::rast` is written for raster. When `return = FALSE`, the file name
#'   or directory (when there're multiple files) is returned.
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
#' dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
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
                          gene_col = "gene", cell_col = "cell_id", z = 3L,
                          phred_col = "qv", min_phred = 20, split_col = NULL,
                          not_in_cell_id = "-1",
                          z_option = c("split", "3d"),
                          file_out = NULL, BPPARAM = SerialParam(),
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
  if (!ext %in% c("csv", "tsv", "txt", "parquet")) {
    stop("The file must be one of csv, tsv, txt, or parquet")
  }
  if (!is.null(file_out)) {
    file_out <- normalizePath(file_out, mustWork = FALSE)
    if (!dir.exists(dirname(file_out)))
      dir.create(dirname(file_out))
    file_dir <- file_path_sans_ext(file_out)
    # File already exists, skip processing
    # make sure if z = "all" transcripts are read from ./detected_transcripts dir
    if (file.exists(file_out) && !dir.exists(file_out) && z != "all") {
        if (!return) return(file_out)
        out <- sfarrow::st_read_parquet(file_out)
        rownames(out) <- out$ID
        return(out)
    } else if (dir.exists(file_dir)) {
      # Multiple files
      pattern <- "\\.parquet$"
      # Need to deal with z-planes
      if (z != "all") {
          pattern <- paste0("_z", paste0(z, collapse = "|"), pattern)
      }
      fns <- list.files(file_dir, pattern, full.names = TRUE)
      if (!length(fns) && length(z) == 1L) {
          pattern <- "\\.parquet$"
          fns <- list.files(file_dir, pattern, full.names = TRUE)
      }
      if (length(fns)) {
          if (!return) return(file_dir)
          out <- lapply(fns, sfarrow::st_read_parquet)
          # add names to a list
          names(out) <- gsub(".parquet", "",
                             x = list.files(file_dir, pattern))
          out <- lapply(out, function(x) {
              # row names are dropped in st_read/write_parquet
              rownames(x) <- x$ID
              return(x)
          })
        return(out)
      }
    }
  }
  if (!is.numeric(z) && z != "all") {
    stop("z must either be numeric or be 'all' indicating all z-planes.")
  }
  if (ext == "parquet") {
    check_installed("arrow")
    mols <- arrow::read_parquet(file) |> data.table::as.data.table()
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
  if (phred_col %in% names(mols)) {
    mols <- mols[mols[[phred_col]] >= min_phred,]
  }
  message(">>> Converting transcript spots to geometry")
  if (dest == "colGeometry") {
      if (!length(cell_col) || any(!cell_col %in% names(mols)))
          stop("Column indicating cell ID not found.")
      mols <- mols[mols[[cell_col[1]]] != not_in_cell_id,]
      if (length(cell_col) > 1L) {
          if (!is.data.table(mols)) ..cell_col <- cell_col
          cell_col_use <- do.call(paste, c(mols[,..cell_col], sep = "_"))
          mols$cell_id_ <- cell_col_use
          mols[,cell_col] <- NULL
          cell_col <- "cell_id_"
      }
  }
  if (z_option == "split" && use_z) {
      mols <- split(mols, mols[[spatialCoordsNames[3]]])
      mols <- lapply(mols, .mols2geo_split, dest = dest,
                     spatialCoordsNames = spatialCoordsNames[1:2],
                     gene_col = gene_col, cell_col = cell_col, BPPARAM = BPPARAM,
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
                              BPPARAM, not_in_cell_id, split_col)
  }

  if (!is.null(file_out)) {
      message(">>> Writing reformatted transcript spots to disk")
      if (is(mols, "sf")) {
          suppressWarnings(sfarrow::st_write_parquet(mols, file_out))
          if (!return) return(file_out)
      } else {
          if (!dir.exists(file_dir)) dir.create(file_dir)
          suppressWarnings({
              bplapply(names(mols), function(n) {
                  name_use <- gsub("/", ".", n)
                  sfarrow::st_write_parquet(mols[[n]],
                                            file.path(file_dir,
                                                      paste0(name_use, ".parquet")))
              }, BPPARAM = SerialParam(progressbar = TRUE))
          })
          if (!return) return(file_dir)
      }
  }
  return(mols)
}

#' @rdname formatTxSpots
#' @export
addTxSpots <- function(sfe, file, sample_id = NULL,
                       spatialCoordsNames = c("global_x", "global_y", "global_z"),
                       gene_col = "gene", z = 3L,
                       phred_col = "qv", min_phred = 20, split_col = NULL,
                       z_option = c("split", "3d"),
                       file_out = NULL, BPPARAM = SerialParam()) {
  sample_id <- .check_sample_id(sfe, sample_id)
  dest <- "rowGeometry"
  mols <- formatTxSpots(file, dest = dest, spatialCoordsNames = spatialCoordsNames,
                        gene_col = gene_col, z = z, phred_col = phred_col,
                        min_phred = min_phred, split_col = split_col,
                        z_option = z_option, file_out = file_out,
                        BPPARAM = BPPARAM)
  if (is(mols, "sf")) {
      txSpots(sfe, withDimnames = TRUE) <- mols
  } else if (is.list(mols)) {
      rowGeometries(sfe) <- c(rowGeometries(sfe), mols)
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
#' by cell compartment. Only relevant when `add_molecules = TRUE`.
#' @return An SFE object. If reading transcript spots (`add_molecules = TRUE`),
#' then the reformatted transcript spots are saved to file specified in the
#' `file_out` argument, which is by default `tx_spots.parquet` in the same
#' directory as the rest of the data.
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
    mat <- fread(fn_mat)
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
    message(">>> Constructing cell polygons")
    polys <- df2sf(polys, spatialCoordsNames = c("x_global_px", "y_global_px"),
                   geometryType = "POLYGON",
                   id_col = "cellID", BPPARAM = BPPARAM)
    polys <- polys[match(meta$cell_ID, polys$ID),]
    sfe <- SpatialFeatureExperiment(list(counts = mat), colData = meta,
                                    spatialCoordsNames = c("CenterX_global_px", "CenterY_global_px"),
                                    unit = "full_res_image_pixel")
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
