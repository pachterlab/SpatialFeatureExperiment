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
  inds <- c(geometries %>% 
              purrr::map(., length) %>% 
              unlist() > 0)
  # remove empty elements
  geometries %<>%
    purrr::keep(., inds) %>%
    purrr::map(., .f = function(m) m %>% t() %>% list() %>% sf::st_polygon())
  # keep non-emplty elements
  df <- st_sf(geometry = sf::st_sfc(geometries),
              ID = cell_ids[inds %>% which()],
              ZIndex = z)
  df
}

#' @importFrom sf st_is_empty
#' @importFrom BiocParallel bplapply
.filter_polygons <- function(polys, min_area, BPPARAM = SerialParam()) {
  # Sanity check on nested polygon lists
  test.segs <- lapply(polys %>% st_geometry() %>% seq, function(i) { 
    polys %>% 
      st_geometry() %>% .[[i]] %>% length() }) %>% unlist()
  if (any(test.segs > 1)) {
    segs.art.index <- which(test.segs > 1)
    warning("Sanity checks on cell segmentation polygons:", "\n", 
            ">>> ..found ", segs.art.index %>% length,
            " cells with (nested) polygon lists", "\n",
            ">>> ..applying filtering")
  }
  # remove empty elements
  polys.orig <- polys
  polys %<>% filter(!st_is_empty(.))
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
                paste(multi_inds %>% head(10), # necessary to print all?
                      collapse = ", "),
                ". The largest piece is kept.")
        which_keep[multi_inds] <- lapply(areas[multi_inds], which.max)
      }
      inds <- lengths(which_keep) > 0L
      polys <- polys[inds,]
      # using parallelization, else can take a while when `which_keep` length is towards 100K
      which_keep <- unlist(which_keep[inds])
      geo <- st_geometry(polys)
      new_geo <- 
        bplapply(seq_along(which_keep), function(i) {
          geo[[i]] <- st_cast(geo[i], "POLYGON")[[which_keep[i]]] %>% 
            unique() %>% # remove any duplicates
            st_polygon()
        }, BPPARAM = BPPARAM) %>% st_sfc()
      st_geometry(polys) <- new_geo
      
    } else if (is.null(min_area)) {
      # use only maximal area, ie the largest polygon
      warning(">>> ..keeping polygons with the largest area only")
      which_keep <- 
        lapply(areas, function(x) which.max(x)) %>% unlist()
      geo <- st_geometry(polys)  
      new_geo <- 
        bplapply(seq_along(which_keep), function(i) {
          geo[[i]] <- st_cast(geo[i], "POLYGON")[[which_keep[i]]] %>% 
            unique() %>% # remove any duplicates
            st_polygon()
        }, BPPARAM = BPPARAM) %>% st_sfc()
      st_geometry(polys) <- new_geo
    }
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
#' @param image Which image(s) to load, can be "DAPI", "PolyT", "Cellbound" or any combination of them.
#' @param min_area Minimum cell area in square microns. Anything smaller will be
#'   considered artifact or debris and removed.
#' @param filter_counts Keep cells with counts \code{> 0}
#' @param add_molecules If to add transcripts coordinates to an object
#' @param use_bboxes If no segmentation output is present, use \code{cell_metadata} to make bounding boxes instead
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
#' 
#' @importFrom sf st_area st_geometry<-
#' @importFrom terra rast ext vect
#' @importFrom BiocParallel bpmapply bplapply
#' @importFrom rlang check_installed
#' @importFrom SpatialExperiment imgData<-
#' @importFrom dplyr
#' @import magrittr
#' 
#' @examples
#' dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
#' sfe <- readVizgen(dir_use, z = 0L, image = "PolyT",
#' flip = "geometry")
#' 
#' ## Filtering of counts, and addition of molecule coordinates..
#' sfe <- readVizgen(dir_use, z = 0L, image = "PolyT", filter_counts = TRUE,
#' add_molecules = TRUE, flip = "geometry")
readVizgen <-
  function(data_dir,
           z = 3L,
           sample_id = "sample01", 
           min_area = 15,
           image = c("DAPI", "PolyT", "Cellbound"),
           flip = c("geometry", "image", "none"),
           max_flip = "50 MB",
           filter_counts = FALSE,
           add_molecules = FALSE,
           use_bboxes = FALSE,
           BPPARAM = SerialParam()) {
    
    # issue message for packages that need to be installed a priori
    pkgs <- c("vroom", "arrow", "sfarrow", "rhdf5",
              "tidyverse", "sf", "BiocParallel", "Matrix")
    lapply(pkgs %>% length %>% seq, function(i) 
    { !requireNamespace(pkgs[i], quietly = TRUE) } ) %>% 
      unlist %>% 
      { if (c(which(.) > 0) %>% any()) 
      { message("Please install ->", "\n",
                paste0("'", pkgs[which(.)], "'", collapse = ", "), " for this function")}
      }
    
    data_dir <- normalizePath(data_dir, mustWork = TRUE)
    flip <- match.arg(flip)
    image <- match.arg(image, several.ok = TRUE)
    rlang::check_installed("vroom")
    if (z < 0 || z > 6) {
      stop("z must be beween 0 and 6 (inclusive).")
    }
    # in some older data, "PolyT" is named "polyT"
    img_fn <-
      file.path(data_dir, "images") %>%
      list.files(ignore.case = TRUE,
                 pattern = paste0(image, collapse = "|")) %>%
      grep(paste0("_z", z), ., value = TRUE)
    img_fn.names <- 
      paste0(c("mosaic_", paste0("_z", z), ".tif"), 
             collapse = "|") %>%
      gsub(., "", img_fn)
    img_fn <- file.path(data_dir, "images", img_fn)
    names(img_fn) <- img_fn.names
    
    do_flip <- .if_flip_img(img_fn, max_flip)
    if_exists <- file.exists(img_fn)
    if (!all(if_exists)) {
      warning("The image file(s) for ", image[!if_exists],
              " don't exist, or have non-standard file name(s).")
      img_fn <- img_fn[if_exists]
    }
    if (!length(img_fn)) flip <- "none"
    else if (!any(do_flip) && flip == "image") flip <- "geometry"
    
    # Use segmentation output from ".parquet" file
    # # check if ".parquet" file is present
    parq <-
      # look in the current directory
      list.files(data_dir,
                 pattern = ".parquet$",
                 full.names = TRUE) %>%
      { if (length(.) == 0) { 
        # look in the sub directory (if nothing is found)
        list.files(data_dir,
                   pattern = ".parquet$",
                   full.names = TRUE,
                   recursive = TRUE)
      } else { (.) }}
    # set to use .parquet" file if present
    use.parquet <- 
      parq %>% length() %>% any
    
    if (use.parquet) {
      message(">>> Cell segmentations are found in `.parquet` file")
      rlang::check_installed("sfarrow")
      if (length(parq) > 1) {
        # eg, if two files are present:
        # `cellpose_micron_space.parquet`
        # `cellpose_mosaic_space.parquet`
        # use µm units
        parq %<>% 
          grep("micron", ., value = TRUE)
      }
      fn <- parq
      # read file and filter to keep selected z section
      polys <- 
        sfarrow::st_read_parquet(fn) %>% 
        filter(., ZIndex == z)
      # filtering cell polygons
      polys <- .filter_polygons(polys, min_area, 
                                BPPARAM = BPPARAM)
      st_geometry(polys) <- "geometry"
      polys$ID <- polys$EntityID
      polys <- polys[,c("ID", "ZIndex", "Type", "ZLevel", "geometry")]
    } else {
      rlang::check_installed("rhdf5")
      fns <- list.files(file.path(data_dir, "cell_boundaries"),
                        "*.hdf5", full.names = TRUE)
      if(length(fns)) {
        message(">>> Reading '.hdf5' files..")
        polys <- bpmapply(.h52poly_fov, fn = fns, SIMPLIFY = FALSE,
                          BPPARAM = BPPARAM,
                          MoreArgs = list(z = z))
        # dplyr is much more efficient than base R rbind
        polys <- if (length(polys) == 1L) polys[[1]] else do.call(bind_rows, polys)
      } else { warning("No '.hdf5' files present, check input directory -> `data_dir`") 
        polys <- NULL
      }
    }
    if (flip == "geometry" && !is.null(polys)) {
      # Flip the coordinates
      mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
      st_geometry(polys) <- st_geometry(polys) * mat_flip
    }
    
    # get count data file
    mat_fn <-
      list.files(data_dir, 
                 pattern = "cell_by_gene",
                 full.names = TRUE) %>%
      { if (length(.) == 0) { 
        list.files(data_dir,
                   pattern = "cell_by_gene",
                   full.names = TRUE,
                   recursive = TRUE)
      } else { (.) }} %>% 
      { if (length(.) > 1) { 
        stop("There are > 1 `cell_by_gene` files",
             "\n", "make sure only 1 file is read")
      } else if (!length(.)) {
        stop("No `cell_by_gene` file is available")
      } else { (.) }}
    # cell ids col is stored in `...1` for older Vizgen data, new data has `cell` col
    suppressMessages(mat <- vroom::vroom(mat_fn, 
                                         col_types = c("c") # convert 1st col to a character
    ))
    
    # get spatial metadata file
    rlang::check_installed("tibble")
    meta_fn <- 
      list.files(data_dir, 
                 pattern = "cell_metadata",
                 full.names = TRUE) %>%
      { if (length(.) == 0) { 
        list.files(data_dir,
                   pattern = "cell_metadata",
                   full.names = TRUE,
                   recursive = TRUE)
      } else { (.) }} %>% 
      { if (length(.) > 1) { 
        stop("There are > 1 `cell_metadata` files",
             "\n", "make sure only 1 file is read")
      } else if (!length(.)) {
        stop("No `cell_metadata` file is available")
      } else { (.) }}
    suppressMessages(metadata <- vroom::vroom(meta_fn, col_types = c("c")))
    metadata %<>% 
      { if (!is.null(polys)) {
        dplyr::slice(., match(polys %>% pull(ID), table = pull(., 1))) 
      } else { (.) }} %>%
      { if ((names(.) == "transcript_count") %>% any() && filter_counts) {
        message(">>> ..filtering `cell_metadata` - keep cells with `transcript_count` > 0")
        filter(., transcript_count > 0)
      } else { (.) }} %>%
      # add cell ids to rownames (this removes cell ids col as well)
      tibble::column_to_rownames(., var = names(.)[1])
    if (flip == "geometry") {
      metadata$center_y <- -metadata$center_y
    }
    
    # convert counts df to sparse matrix
    m <-
      mat %>% 
      { if (!is.null(polys)) {
        dplyr::slice(., match(polys %>% pull(ID), table = pull(., 1))) 
      } else { (.) }} %>%
      tibble::column_to_rownames(., var = names(.)[1]) %>%
      # add filtering or match to (filtered) metadata
      { if ((names(metadata) == "transcript_count") %>% any() 
            && filter_counts) {
        # match cell ids to filtered metadata
        filter(., rownames(.) %in% rownames(metadata))
      } else if (filter_counts) { # filter only count matrix
        message(">>> filtering `cell_by_gene` - keep cells with counts > 0")
        filter_all(., any_vars(. > 0))
      } else { # pass it further as it is
        (.) }
      } %>%
      as.matrix() %>%
      as("CsparseMatrix") %>% # convert to sparse matrix
      Matrix::t() # transpose sparse matrix
    
    # if only count matrix is filtered..
    #..match cell IDs from filtered count matrix to cell metadata df
    if (!(names(metadata) == "transcript_count") %>% any() && filter_counts) {
      metadata %<>% filter(rownames(.) %in% colnames(m))
    }
    
    # check matching cell ids in geometries
    if (!is.null(polys) &&
        !identical(polys %>% pull(ID), 
                   m %>% colnames)) {
      # filter geometries as well
      matched.cells <- match(m %>% colnames, polys %>% pull(ID))
      message(">>> filtering geometries to match ", length(matched.cells), 
              " cells with counts > 0") 
      polys %<>% dplyr::slice(matched.cells)
    }
    
    if (any(if_exists)) {
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
    }
    
    # add molecule coordiantes to rowGeometries
    if (add_molecules) {
      rlang::check_installed("purrr")
      
      # get molecule coordiantes file 
      mols_fn <- 
        list.files(data_dir, 
                   pattern = "detected_transcripts",
                   full.names = TRUE) %>%
        { if (length(.) == 0) { 
          list.files(data_dir,
                     pattern = "detected_transcripts",
                     full.names = TRUE,
                     recursive = TRUE)
        } else { (.) }} %>% 
        { if (length(.) > 1) { 
          stop("There are > 1 `detected_transcripts` files",
               "\n", "make sure only 1 file is read")
        } else if (!length(.)) {
          stop("No `detected_transcripts` file is available")
        } else { (.) }}
      suppressMessages(mols <- vroom::vroom(mols_fn, 
                                            col_types = vroom::cols(cell_id = "c")))
      # filtering..
      mols %<>%
        # match cell id (if present) with count matrix
        { if ((names(.) == "cell_id") %>% any) {
          # keep only selected z-plane
          # remove transcripts that are not associated with a cell
          filter(., global_z == z, !cell_id == "-1") %>%
            dplyr::slice(pull(., cell_id) %>% match(., m %>% colnames))
        } else { 
          # keep only selected z-plane
          filter(., global_z == z) }} %>%
        # keep x,y coords (in µm), gene and cell_id cols
        select(contains(c("global_x", "global_y",
                          "gene", "cell_id"))) %>%
        # rename to standard x and y
        rename(x = global_x, y = global_y) %>%
        as.data.frame() %>%
        # split by gene
        group_by(gene) %>%
        group_split() %>% 
        purrr::map(., .f = function(i) st_sf(gene = i %>% pull(gene) %>% unique,
                                             geometry =              
                                               i %>% 
                                               select(x, y) %>%
                                               as.matrix() %>%
                                               st_multipoint() %>%
                                               st_sfc() %>%
                                               st_cast('MULTIPOINT'))) %>%
        do.call("bind_rows", .)
      
      # add gene names to rows
      rownames(mols) <- mols %>% pull(gene)
      # make sure genes names correspond to rows of count matrix, ie same order
      mols %<>%
        dplyr::slice(., match(m %>% rownames(), table = rownames(.)))
    }
    
    # Takes a while to make the POINT geometry for the centroids, not too bad
    sfe <- SpatialFeatureExperiment(assays = list(counts = m),
                                    colData = metadata,
                                    sample_id = sample_id,
                                    spatialCoordsNames = c("center_x", "center_y"),
                                    unit = "micron")
    
    # If none of segmentations are present, make bounding boxes
    # NOTE: might take some time to run
    if (use_bboxes && is.null(polys)) {
      message(">>> Creating bounding boxes from `cell_metadata`")
      bboxes <-
        bplapply(metadata %>% nrow() %>% seq(), 
                 function(cell.id) {
                   st_sf(ID = sfe %>% colnames() %>% .[cell.id],
                         geometry = select(slice(metadata, cell.id), 
                                           starts_with(c("min_", "max_"))) %>%
                           setNames(c("xmin", "ymin", 
                                      "xmax", "ymax")) %>%
                           unlist() %>%
                           st_bbox() %>% st_as_sfc())
                 }, BPPARAM = BPPARAM) %>% 
        do.call("bind_rows", .)
      
      # add gene names to rows
      rownames(bboxes) <- bboxes %>% pull(ID)
      # make sure genes names correspond to rows of count matrix, ie same order
      bboxes %<>%
        dplyr::slice(., match(sfe %>% colnames(), table = rownames(.)))
      bboxes$ID <- NULL
      cellSeg(sfe) <- bboxes
    }
    
    if (!is.null(polys)) {
      rownames(polys) <- polys$ID
      polys$ID <- NULL
      cellSeg(sfe) <- polys
    }
    
    if (any(if_exists)) { imgData(sfe) <- img_df }
    
    if (add_molecules) { 
      # store molecules in rowGeometries
      message(">>> Storing molecule coordinates in `rowGeometries`") 
      suppressWarnings(rowGeometries(sfe, withDimnames = FALSE) <- mols %>% list())
      dimGeometryNames(sfe, MARGIN = 1) <- "molecules"
    }
    sfe
  }
