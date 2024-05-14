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

.make_sql_query <- function(fn, gene_select, gene_col) {
    if (is.null(gene_select)) return(NA)
    lyr_name <- basename(fn) |> file_path_sans_ext()
    gene_part <- paste0("('", paste0(gene_select, collapse = "','"), "')")
    paste0("SELECT * FROM ", lyr_name, " WHERE ", gene_col, " IN ", gene_part)
}

.read_tx_output <- function(file_out, z, z_option, gene_col, return) {
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
            out <- lapply(fns, sfarrow::st_read_parquet)
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
        out <- sfarrow::st_read_parquet(file_out)
        rownames(out) <- out[[gene_col]]
        return(out)
    }
}

#' Read transcript spots of select genes
#'
#' I speculate that in practice, the most common use of the transcript spots is
#' visualization, and only a few genes can be visualized at a time or the spots
#' will overcrowd. Then it doesn't make sense to load the transcript spots of
#' all genes into memory as they can take up a lot of memory. The function
#' \code{readSelectTx} reads transcript spots of select genes into R, and the
#' function \code{addSelectTx} adds them to \code{rowGeometries} of the SFE
#' object.
#'
#' @inheritParams formatTxSpots
#' @param sfe A `SpatialFeatureExperiment` object.
#' @param file File path of a GeoParquet file (e.g. already reformatted with the
#'   \code{\link{formatTxSpots}} or \code{\link{addTxSpots}} function, should
#'   have already flipped to match image if necessary).
#' @param gene_select Character vector of a subset of genes. If \code{NULL},
#'   then all genes that have transcript spots are added. Only relevant when
#'   reading data from formatted files on disk. If specified, then \code{return
#'   = TRUE}.
#' @param swap_rownames Name of a column in \code{rowData(sfe)} to use as gene
#'   identifiers in place of the actual row names. In some cases this may be
#'   needed to match each transcript spot MULTIPOINT geometry to rows of
#'   \code{sfe}.
#' @note The GDAL Parquet driver is required for this function, though not for
#'   other functions that work with GeoParquet files. GDAL Parquet driver has
#'   been supported since GDAL 3.5.0, but is not part of the default
#'   installation. The \code{z} and \code{z_option} arguments are there since
#'   the file names contain z-plane information when relevant.
#' See the \href{https://gdal.org/drivers/vector/parquet.html}{GDAL documentation
#' page for the Parquet driver}.
#' @return When there are multipel parquet files to be read, a list of sf data
#'   frames with MULTIPOINT geometry for genes selected. When there is only one
#'   file, then one sf data frame. For \code{addSelectTx}, an SFE object with
#'   the transcript spots of the selected genes added.
#' @name readSelectTx
#' @concept Transcript spots
#' @export
#' @examples
#' library(SFEData)
#' if (gdalParquetAvailable()) {
#'     fp <- tempdir()
#'     dir_use <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
#'     fn_tx <- formatTxTech(dir_use, tech = "Xenium", flip = TRUE, return = FALSE,
#'                           file_out = file.path(dir_use, "tx_spots.parquet"))
#'     gene_select <- c("ACE2", "BMX")
#'     df <- readSelectTx(fn_tx, gene_select)
#'
#'     sfe <- readXenium(dir_use)
#'     sfe <- addSelectTx(sfe, fn_tx, head(rownames(sfe), 5), swap_rownames = "Symbol")
#'     unlink(dir_use, recursive = TRUE)
#' }
readSelectTx <- function(file, gene_select, z = "all",
                         z_option = c("3d", "split")) {
    if (!gdalParquetAvailable()) {
        stop("GDAL Parquet driver is required to selectively read genes.")
    }
    file <- normalizePath(file, mustWork = FALSE)
    file_dir <- file_path_sans_ext(file)
    z_option <- match.arg(z_option)
    gene_col <- "gene"
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
            out <- lapply(fns, function(x) {
                q <- .make_sql_query(x, gene_select, gene_col)
                out <- st_read(x, query = q, int64_as_string = TRUE, quiet = TRUE,
                               crs = NA)
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
    } else if (file.exists(file) && !dir.exists(file_dir)) {
        q <- .make_sql_query(file, gene_select, gene_col)
        out <- st_read(file, query = q, int64_as_string = TRUE, quiet = TRUE,
                       crs = NA)
        rownames(out) <- out[[gene_col]]
        return(out)
    }
}

#' @rdname readSelectTx
#' @export
addSelectTx <- function(sfe, file, gene_select, sample_id = 1L,
                        z = "all", z_option = c("3d", "split"),
                        swap_rownames = NULL) {
    sample_id <- .check_sample_id(sfe, sample_id)
    z_option <- match.arg(z_option)
    gene_select <- .id2symbol(sfe, gene_select, swap_rownames)
    mols <- readSelectTx(file, gene_select, z, z_option)
    rownames(mols) <- .symbol2id(sfe, rownames(mols), swap_rownames)
    if (!is(mols, "sf")) {
        rowGeometries(sfe, sample_id = sample_id, partial = TRUE) <- mols
    } else {
        txSpots(sfe, sample_id, partial = TRUE) <- mols
    }
    sfe
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
#' @return A sf data frame for vector geometries if `file_out` is not set.
#'   `SpatRaster` for raster. If there are multiple files written, such as when
#'   splitting by cell compartment or when `dest = "colGeometry"`, then a
#'   directory with the same name as `file_out` will be created (but without the
#'   extension) and the files are written to that directory with informative
#'   names. `parquet` files that can be read with `st_read` is written for
#'   vector geometries. When `return = FALSE`, the file name or directory
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
#' @concept Transcript spots
#' @return The `sf` data frame, or path to file where geometries are written if
#'   `return = FALSE`.
#' @rdname formatTxSpots
#' @examples
#' # Default arguments are for MERFISH
#' fp <- tempdir()
#' dir_use <- SFEData::VizgenOutput(file_path = file.path(fp, "vizgen_test"))
#' g <- formatTxSpots(file.path(dir_use, "detected_transcripts.csv"))
#' unlink(dir_use, recursive = TRUE)
#'
#' # For CosMX, note the colnames, also dest = "colGeometry"
#' # Results are written to the tx_spots directory
#' dir_use <- SFEData::CosMXOutput(file_path = file.path(fp, "cosmx_test"))
#' cg <- formatTxSpots(file.path(dir_use, "Run5642_S3_Quarter_tx_file.csv"),
#' dest = "colGeometry", z = "all",
#' cell_col = c("cell_ID", "fov"),
#' gene_col = "target", not_in_cell_id = "0",
#' spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
#' file_out = file.path(dir_use, "tx_spots"))
#' # Cleanup
#' unlink(dir_use, recursive = TRUE)
formatTxSpots <- function(file, dest = c("rowGeometry", "colGeometry"),
                          spatialCoordsNames = c("global_x", "global_y", "global_z"),
                          gene_col = "gene", cell_col = "cell_id", z = "all",
                          phred_col = "qv", min_phred = 20, split_col = NULL,
                          not_in_cell_id = c("-1", "UNASSIGNED"),
                          z_option = c("3d", "split"), flip = FALSE,
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
    if (!ext %in% c("csv", "gz", "tsv", "txt", "parquet")) {
        stop("The file must be one of csv, gz, tsv, txt, or parquet")
    }
    if (!is.null(file_out) && (file.exists(file_out) ||
                               dir.exists(file_path_sans_ext(file_out)))) {
        out <- .read_tx_output(file_out, z, z_option, "gene", return)
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
    names(mols)[names(mols) == gene_col] <- "gene"
    gene_col <- "gene"
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
            suppressWarnings(sfarrow::st_write_parquet(mols, file_out))
            if (!return) return(file_out)
        } else {
            if (!dir.exists(file_dir)) dir.create(file_dir)
            suppressWarnings({
                bplapply(names(mols), function(n) {
                    name_use <- gsub("/", ".", n)
                    suppressWarnings(sfarrow::st_write_parquet(mols[[n]], file.path(file_dir, paste0(name_use, ".parquet"))))
                }, BPPARAM = SerialParam(progressbar = TRUE))
            })
            if (!return) return(file_dir)
        }
    }
    return(mols)
}

#' @rdname formatTxSpots
#' @export
addTxSpots <- function(sfe, file, sample_id = 1L,
                       spatialCoordsNames = c("global_x", "global_y", "global_z"),
                       gene_col = "gene", z = "all",
                       phred_col = "qv", min_phred = 20, split_col = NULL,
                       z_option = c("3d", "split"), flip = FALSE,
                       file_out = NULL, BPPARAM = SerialParam()) {
    sample_id <- .check_sample_id(sfe, sample_id)
    z_option <- match.arg(z_option)
    dest <- "rowGeometry"
    mols <- formatTxSpots(file, dest = dest, spatialCoordsNames = spatialCoordsNames,
                          gene_col = gene_col, z = z, phred_col = phred_col,
                          min_phred = min_phred, split_col = split_col,
                          flip = flip, z_option = z_option, file_out = file_out,
                          BPPARAM = BPPARAM, return = TRUE)
    if (is(mols, "sf")) {
        txSpots(sfe, withDimnames = TRUE) <- mols
    } else if (is(mols, "list")) {
        rowGeometries(sfe) <- mols
    }

    # make sure that sfe and rowGeometries have the same features
    # NOTE, if `min_phred = NULL`, no filtering of features occur
    if (!is.null(min_phred)) {
        gene_col <- "gene" # It's always "gene" from formatTxSpots
        gene_names <-
            lapply(rowGeometries(sfe), function(i) {
                na.omit(i[[gene_col]])
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
    }
    sfe
}

#' Read and process transcript spots for specific commercial technologies
#'
#' To preset parameters such as \code{spatialCoordsNames}, \code{gene_col},
#' \code{cell_col}, and \code{phred_col} that are standard for the output of the
#' technology.
#'
#' @inheritParams formatTxSpots
#' @inheritParams readVizgen
#' @param tech Which technology whose output to read, must be one of "Vizgen",
#' "Xenium", or "CosMX" though more technologies may be added later.
#' @param z Which z-planes to read. Always "all" for Xenium where the z
#'   coordinates are not discrete.
#' @param split_cell_comps Only relevant to CosMX whose transcript spot file
#' assigns the spots to cell components. Setting this argument to \code{TRUE}
#' @return The `sf` data frame, or path to file where geometries are written if
#'   `return = FALSE`.
#' @name formatTxTech
#' @concept Transcript spots
#' @export
#' @examples
#' library(SFEData)
#' fp <- tempdir()
#' dir_use <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
#' fn_tx <- formatTxTech(dir_use, tech = "Xenium", flip = TRUE, return = FALSE,
#'                       file_out = file.path(dir_use, "tx_spots.parquet"))
#'
formatTxTech <- function(data_dir, tech = c("Vizgen", "Xenium", "CosMX"),
                         dest = c("rowGeometry", "colGeometry"),
                         z = "all", min_phred = 20,
                         split_cell_comps = FALSE,
                         z_option = c("3d", "split"), flip = FALSE,
                         file_out = NULL, BPPARAM = SerialParam(),
                         return = TRUE) {
    data_dir <- normalizePath(data_dir, mustWork = TRUE)
    tech <- match.arg(tech)
    c(spatialCoordsNames, gene_col, cell_col, fn) %<-%
        .get_tech_tx_fields(tech, data_dir)
    if (tech == "CosMX" && split_cell_comps)
        split_col <- "CellComp"
    else split_col <- NULL
    # TODO: check if z is valid here for all technologies. Need new internal function.
    if (tech == "Xenium") {
        z <- "all"
        z_option <- "3d"
    }
    formatTxSpots(fn, dest = dest, spatialCoordsNames = spatialCoordsNames,
                  gene_col = gene_col, cell_col = cell_col, z = z,
                  min_phred = min_phred, split_col = split_col,
                  z_option = z_option, flip = flip,
                  file_out = file_out, BPPARAM = BPPARAM,
                  return = return
    )
}

#' @rdname formatTxTech
#' @export
addTxTech <- function(sfe, data_dir, sample_id = 1L,
                      tech = c("Vizgen", "Xenium", "CosMX"),
                      z = "all", min_phred = 20,
                      split_cell_comps = FALSE,
                      z_option = c("3d", "split"), flip = FALSE,
                      file_out = NULL, BPPARAM = SerialParam()) {
    c(spatialCoordsNames, gene_col, cell_col, fn) %<-%
        .get_tech_tx_fields(tech, data_dir)
    if (tech == "CosMX") flip <- FALSE
    if (tech == "CosMX" && split_cell_comps)
        split_col <- "CellComp"
    else split_col <- NULL
    if (tech == "Xenium") {
        z <- "all"
        z_option <- "3d"
    } else min_phred <- NULL # So far only Xenium has phred score
    addTxSpots(sfe, file = fn, sample_id = sample_id,
               spatialCoordsNames = spatialCoordsNames,
               gene_col = gene_col, z = z,
               phred_col = "qv", min_phred = min_phred, split_col = split_col,
               z_option = z_option, flip = flip,
               file_out = file_out, BPPARAM = BPPARAM)
}
