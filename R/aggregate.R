#' Aggregate transcript spots from file
#'
#' This function reads the transcript spot file from the standard output of the
#' commercial technologies (not GeoParquet) for spatial aggregation where the
#' spots are assigned to polygons such as cells or spatial bins. Presets for
#' Xenium, MERFISH, and CosMX are available. For Vizgen and Xenium, the images
#' can be added when \code{add_images = TRUE}.
#'
#' @param df If the file is already loaded into memory, a data frame (sf) with
#'   columns for the x, y, and optionally z coordinates and gene assignment of
#'   each transcript spot. If specified, then argument \code{file} will be
#'   ignored.
#' @param by A \code{sfc} or \code{sf} object for spatial aggregation.
#' @param new_geometry_name Name to give to the new \code{colGeometry} in the
#'   output. Defaults to "bins".
#' @param flip_geometry Logical, whether to flip the transcript spot geometries
#'   to match the images if added later.
#' @param image String, which image(s) to add to the output SFE object. Not
#'   applicable to CosMX. See \code{\link{readVizgen}} and
#'   \code{\link{readXenium}} for options and multiple images can be specified.
#'   If \code{NULL}, then the default from the read function for the technology
#'   will be used.
#' @param sparse Logical, whether the gene count matrix from aggregating
#'   transcript spots should be sparse. When the bins are large, the matrix will
#'   not be very sparse so using sparse matrix will not save memory, but when
#'   the bins are small, sparsity is worth it.
#' @param BPPARAM bpparam object to specify parallel computing over genes. If a
#'   lot of memory is used, then stick to `SerialParam()`.
#' @inheritParams formatTxTech
#' @inheritParams formatTxSpots
#' @inheritParams readXenium
#' @inheritParams sf::st_make_grid
#' @inheritParams SpatialFeatureExperiment
#' @note The resulting SFE object often includes geometries (e.g. grid cells)
#'   outside tissue, because there can be transcript spots detected outside the
#'   tissue. Also, bins at the edge of the tissue that don't fully overlap with
#'   the tissue will have lower transcript counts; this may have implications to
#'   downstream spatial analyses.
#' @return A SFE object with count matrix for number of spots of each gene in
#'   each geometry. Geometries with no spot are removed.
#' @importFrom data.table as.data.table
#' @importFrom sf st_make_grid
#' @concept Geometric operations
#' @export
aggregateTx <- function(file, df = NULL, by = NULL, sample_id = "sample01",
                        spatialCoordsNames = c("X", "Y", "Z"),
                        gene_col = "gene",
                        phred_col = "qv", min_phred = 20, flip_geometry = FALSE,
                        cellsize = NULL, square = TRUE, flat_topped = FALSE,
                        new_geometry_name = "bins", unit = "micron", sparse = FALSE,
                        BPPARAM = SerialParam()) {
    # This is only for one file, one sample
    if (!is.null(df)) file <- df
    mols <- .check_tx_file(file, spatialCoordsNames, gene_col, phred_col,
                           min_phred, flip_geometry)
    if (inherits(mols, "data.table"))
        mols <- mols[,c(spatialCoordsNames, "gene"), with=FALSE]
    else mols <- mols[,c(spatialCoordsNames, "gene")]
    mols <- df2sf(mols, spatialCoordsNames = spatialCoordsNames,
                  geometryType = "POINT")
    if (is.null(by))
        by <- st_make_grid(mols, cellsize = cellsize, square = square,
                           flat_topped = flat_topped)
    else if (inherits(by, "sf")) by <- st_geometry(by)
    mols <- split(st_geometry(mols), mols[["gene"]])
    if (sparse) { # TODO: try Arrow dataset, querying one gene at a time, 
        # then create TileDB right from the beginning
        # Iterate over the genes, count number of transcripts in each bin for the gene
        ml <- bplapply(seq_along(mols), function(i) {
            x <- mols[[i]]
            ll <- lengths(st_intersects(by, x))
            j <- which(ll > 0)
            data.frame(i = i, j = j, x = ll[j])
        }, BPPARAM = BPPARAM)
        ml <- data.table::rbindlist(ml)
        new_mat <- sparseMatrix(i = ml$i, j = ml$j, x = ml$x, dims = c(length(mols), length(by)),
                                dimnames = list(names(mols), seq_along(by)))
    } else {
        ml <- bplapply(mols, function(x) {
            inds <- st_intersects(by, x)
            lengths(inds)
        }, BPPARAM = BPPARAM)
        new_mat <- matrix(unlist(ml), nrow = length(by), ncol = length(mols),
                          dimnames = list(seq_along(by), names(mols)))
        new_mat <- t(new_mat)
    }
    new_mat <- new_mat[,colSums(new_mat) > 0] # Remove empty grid cells
    grid_sf <- st_sf(geometry = by)
    cgs <- list(bins = grid_sf[colnames(new_mat), "geometry"])
    names(cgs) <- new_geometry_name
    SpatialFeatureExperiment(assays = list(counts = new_mat),
                             colGeometries = cgs, unit = unit)
}

#' @rdname aggregateTx
#' @export
aggregateTxTech <- function(data_dir, df = NULL, by = NULL,
                            tech = c("Vizgen", "Xenium", "CosMX"),
                            sample_id = "sample01",
                            image = NULL,
                            min_phred = 20, flip = c("geometry", "image", "none"),
                            max_flip = "50 MB",
                            cellsize = NULL, square = TRUE, flat_topped = FALSE,
                            new_geometry_name = "bins", sparse = FALSE,
                            BPPARAM = SerialParam()) {
    tech <- match.arg(tech)
    flip <- match.arg(flip)
    c(spatialCoordsNames, gene_col, cell_col, fn) %<-%
        getTechTxFields(tech, data_dir)
    if (tech == "Xenium") {
        c(xoa_version, major_version, minor_version, instrument_version) %<-%
            .get_XOA_version(data_dir)
    }
    img_choices <- image <- switch (
        tech,
        Vizgen = c("DAPI", "PolyT", "Cellbound"),
        Xenium = if (major_version > 1L) "morphology_focus" else c("morphology_focus", "morphology_mip"),
        NA
    )
    if (is.null(image)) {
        image <- img_choices
    } else {
        image <- match.arg(image, img_choices, several.ok = TRUE)
    }
    if (tech == "Xenium") {
        c(img_df, flip) %<-% .get_xenium_images(data_dir, image, major_version,
                                                flip, max_flip, sample_id)
    } else if (tech == "Vizgen") {
        c(img_df, flip) %<-% .get_vizgen_images(data_dir, image, flip, max_flip,
                                                z = "all", sample_id = sample_id)
    } else {
        img_df <- NULL
        flip <- "none"
    }
    sfe <- aggregateTx(file = fn, df = df, by = by, sample_id = sample_id,
                spatialCoordsNames = spatialCoordsNames,
                gene_col = gene_col,
                phred_col = "qv", min_phred = min_phred,
                flip_geometry = (flip == "geometry"),
                cellsize = cellsize, square = square, flat_topped = flat_topped,
                new_geometry_name = new_geometry_name)
    imgData(sfe) <- img_df
    sfe
}

.aggregate_num <- function(mat, inds, FUN, fun_name, BPPARAM) {
    if (fun_name %in% c("sum", "mean")) {
        i <- unlist(inds) # cell index
        ll <- lengths(inds) # bin index
        j <- rep(seq_along(inds), times = ll)
        if (fun_name == "sum")
            mat_agg <- sparseMatrix(i = i, j = j, x = 1,
                                    dims = c(ncol(mat), length(inds)))
        else if (fun_name == "mean") {
            mat_agg <- sparseMatrix(i = i, j = j,
                                    x = rep(1/ll, times = ll),
                                    dims = c(ncol(mat), length(inds)))
        }
        # mat has genes in rows and cells in columns
        out_agg <- mat %*% mat_agg
        # Variance and standard deviation can also be expressed with matrix
        # algebra. Not implementing yet since most likely they're not commonly
        # used for aggregation purposes.
    } else if (fun_name %in% getNamespaceExports("sparseMatrixStats")) {
        out_agg <- bplapply(inds, function(i) {
            m <- mat[,i,drop = FALSE]
            FUN(m)
        }, BPPARAM = BPPARAM)
        out_agg <- matrix(unlist(out_agg), ncol = length(inds))
        rownames(out_agg) <- rownames(mat)
    } else stop("Function ", fun_name, " not supported for aggregating SFE.")
    if (inherits(out_agg, "dgeMatrix")) out_agg <- as.matrix(out_agg)
    out_agg
}

# For one sample
.aggregate_sample_cell <- function(x, by, FUN, fun_name, colGeometryName, join,
                                   new_geometry_name, BPPARAM) {
    cg <- colGeometry(x, type = colGeometryName)
    # Preliminary check of overlap
    if (!join(st_as_sfc(st_bbox(by)), st_as_sfc(st_bbox(cg)), sparse = FALSE))
        stop("`by` does not overlap with this sample")
    inds <- join(by, cg) # somewhat faster than join(cg, by) when by has fewer geometries than cg
    not_empty <- which(lengths(inds) > 0L)
    if (!length(not_empty)) stop("`by` does not overlap with this sample")
    by <- by[not_empty]
    inds <- inds[not_empty]
    cnts <- counts(x)
    # Aggregate the gene count matrix, only do the counts assay
    # Express this as a matrix multiplication for sum. Use sweep to deal with mean
    cnts_agg <- .aggregate_num(cnts, inds, FUN, fun_name, BPPARAM)
    # Aggregate numeric columns of colData; logical are converted to numeric
    df <- colData(x)
    df$sample_id <- NULL
    which_logical <- which(vapply(df, is.logical, FUN.VALUE = logical(1)))
    for (i in which_logical) {
        df[,i] <- as.integer(df[,i])
    }
    numeric_cols <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    if (any(numeric_cols)) {
        cd <- t(as.matrix(df[,numeric_cols, drop = FALSE]))
        cd_agg <- .aggregate_num(cd, inds, FUN, fun_name, BPPARAM) |> Matrix::t()
        cd_agg <- as(cd_agg, "DFrame")
    } else cd_agg <- make_zero_col_DFrame(length(inds))
    # Deal with anything that is neither numerical nor logical
    # Primarily character and factor. What about list columns? I don't know.
    # For now, the list columns will be split just like categorical variables
    # and become nested lists
    if (!all(numeric_cols)) {
        df_bin <- data.frame(bin = rep(seq_along(inds), times = lengths(inds)),
                             index = unlist(inds))
        df_bin <- df_bin[order(df_bin$index),]
        names_not_num <- names(df)[!numeric_cols]
        cat_agg <- matrix(NA, nrow = length(inds), ncol = length(names_not_num))
        colnames(cat_agg) <- names_not_num
        cat_agg <- data.frame(cat_agg)
        if (nrow(df_bin) != nrow(df)) {
            df_inds <- data.frame(index = seq_len(nrow(df)))
            df_bin <- merge(df_inds, df_bin, all.x = TRUE, by = "index")
        }
        for (n in names_not_num)
            cat_agg[[n]] <- split(df[[n]], df_bin$bin)
        cd_agg <- cbind(cat_agg, cd_agg)
    }
    cgs <- list(bins = st_sf(geometry = by))
    names(cgs) <- new_geometry_name
    out <- SpatialFeatureExperiment(assays = list(counts = cnts_agg), colData = cd_agg,
                                    sample_id = sampleIDs(x),
                                    colGeometries = cgs)
    colnames(out) <- seq_along(by)
    rownames(out) <- rownames(x)
    out
}

# Might turn this into an exported function
.aggregate_sample_tx <- function(x, by, rowGeometryName, new_geometry_name, 
                                 sparse = FALSE,
                                 BPPARAM = SerialParam()) {
    rg <- rowGeometry(x, rowGeometryName)
    if (!is.null(st_z_range(rg)))
        by <- st_zm(by, drop = FALSE, what = "Z")
    by <- by[lengths(st_intersects(by, rg)) > 0]
    grid_sf <- st_sf(grid_id = seq_along(by), geometry = by)
    tx_coords <- st_coordinates(rg) |> as.data.frame()
    if (ncol(tx_coords) > 3L) scn <- c("X", "Y", "Z") else scn <- c("X", "Y")
    out <- aggregateTx(df = tx_coords, spatialCoordsNames = scn,
                       gene_col = "L1", by = by, sparse = sparse,
                       BPPARAM = BPPARAM)
    rownames(out) <- rownames(x)
    out
}

.aggregate_sample <- function(x, by = NULL, FUN = sum, fun_name,
                              colGeometryName = 1L, rowGeometryName = NULL,
                              join = st_intersects, new_geometry_name = "bins",
                              sparse = FALSE,
                              BPPARAM = SerialParam()) {
    # Here x is an SFE object with one sample
    # by is sfc
    # Can't do S4 method with signature for `by` because the argument `by` isn't
    # in the generic and I don't want to mess with the `aggregate` function in
    # other packages
    if (!is.null(rowGeometryName)) {
        .aggregate_sample_tx(x, by, rowGeometryName, new_geometry_name, 
                             sparse = sparse, BPPARAM = BPPARAM)
    } else {
        .aggregate_sample_cell(x, by, FUN, fun_name, colGeometryName, join,
                               new_geometry_name, BPPARAM)
    }
}

.aggregate_SFE <- function(x, by = NULL, FUN = sum, sample_id = "all",
                           colGeometryName = 1L, rowGeometryName = NULL,
                           cellsize = NULL, square = TRUE, flat_topped = FALSE,
                           new_geometry_name = "bins", join = st_intersects,
                           sparse = FALSE,
                           BPPARAM = SerialParam()) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    if (is.null(by) && is.null(cellsize)) {
        stop("Either `by` or `cellsize` must be specified.")
    }
    # Make grid for multiple samples if `by` is not specified
    if (is.null(by)) {
        by <- .make_grid_samples(x, sample_id,
                                 cellsize, square, flat_topped)
    }
    if (is.list(by) && !inherits(by, "sfc") && !inherits(by, "sf")) {
        if (!any(sample_id %in% names(by)))
            stop("None of the geometries in `by` correspond to sample_id")
        by <- by[intersect(sample_id, names(by))]
    } else {
        if (!inherits(by, "sfc") && !inherits(by, "sf"))
            stop("`by` must be either sf or sfc.")
        if (length(sample_id) > 1L) {
            if (inherits(by, "sfc") || !"sample_id" %in% names(by))
                stop("`by` must be an sf data frame with a column `sample_id`")
            by <- split(st_geometry(by), by$sample_id)
        } else if (inherits(by, "sf")) {
            by <- st_geometry(by)
        }
    }
    if (inherits(by, "sfc")) by <- setNames(list(by), sample_id)
    fun_name <- as.character(substitute(FUN))
    sfes <- splitSamples(x) # Output list should have sample IDs as names
    sfes <- lapply(sample_id, function(s) {
        .aggregate_sample(sfes[[s]], by = by[[s]], FUN = FUN,
                          colGeometryName = colGeometryName,
                          rowGeometryName = rowGeometryName,
                          join = join, fun_name = fun_name,
                          new_geometry_name = new_geometry_name,
                          sparse = sparse, BPPARAM = BPPARAM)
    })
    out <- do.call(cbind, sfes)
    # Add the original rowGeometries back
    rowGeometries(out) <- rowGeometries(x, sample_id = sample_id)
    # Keep imgData
    id_orig <- imgData(x)
    imgData(out) <- id_orig[id_orig$sample_id %in% sample_id,]
    out
}

#' Aggregate data in SFE using geometry
#'
#' Gene expression and numeric columns of \code{colData} will be aggregated with
#' the function specified in \code{FUN}, according to another geometry supplied
#' and a geometry predicate (such as \code{st_intersects}). For example, when
#' the predicate is \code{st_intersects} and a spatial grid is used to
#' aggregate, then the data associated with all cells that intersect with each
#' grid cell will be aggregated with \code{FUN}, such as \code{mean} or
#' \code{sum}. The categorical columns will be collected into list columns, and
#' logical columns will be converted into numeric before applying \code{FUN}.
#'
#' For smFISH-based data where the transcript spots are available, the
#' transcript spots can be used instead of cells to aggregate the gene count
#' matrix, in which case all assays other than \code{counts} will be dropped and
#' \code{FUN} only applies to \code{colData} because the transcript spots are
#' simply counted.
#'
#' What this function does is similar to
#' \href{https://github.com/JEFworks-Lab/SEraster}{SEraster} but more general
#' because any geometry and more aggregation function can be used, not just
#' regular grids, and the aggregation can be performed on the transcript spots.
#'
#' @inheritParams sf::st_make_grid
#' @inheritParams sf::aggregate
#' @param x An SFE object to be aggregated.
#' @param by A \code{sf} data frame whose geometry column is used for
#'   aggregation or \code{sfc} or for multiple samples a list of \code{sfc}
#'   whose names are the sample IDs. For multiple samples, the \code{sf} data
#'   frame must have a column \code{sample_id} to indicate which geometry for
#'   which sample. This argument is optional if \code{cellsize} is specified.
#' @param FUN Function to aggregate the numerical columns in \code{colData} and
#'   the gene count matrix. This can be \code{sum}, \code{mean}, or any function
#'   that takes a numeric matrix as input and returns a numeric vector whose
#'   length is same as the number of rows in the input matrix, such as
#'   \code{rowMedians}. See package \code{matrixStats}. Depending on the
#'   function used for aggregation, numeric columns of \code{colData} may need
#'   to be interpreted differently after aggregation. Aggregation is not done
#'   when aggregating by transcript spots in \code{rowGeometry}. When it's sum
#'   or mean, matrix multiplication is used for aggregation rather than calling
#'   the sum or mean function itself; this is much faster than looping through
#'   the bins and calling the function on each of them.
#' @param sample_id Which samples to aggregate, defaults to "all".
#' @param colGeometryName Which \code{colGeometry} to spatially aggregate the
#'   data, by default the first one.
#' @param rowGeometryName Which \code{rowGeometry} to spatially aggregate
#' @param new_geometry_name Name to give to the new \code{colGeometry} in the
#'   output. Defaults to "bins".
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying parallel
#'   computing when aggregating data with functions other than sum and mean when
#'   aggregating cells. When aggregating transcript spots, this specifies
#'   parallel computing over genes. Defaults to \code{SerialParam()}.
#' @return An SFE object with \code{colGeometry} the same as the geometry
#'   specified in \code{by} or same as the grid specified in \code{cellsize}.
#'   \code{rowGeometries} and \code{rowData} remain the same as in the input
#'   \code{x}. \code{reducedDims}, \code{localResults}, \code{colFeatureData}
#'   (and its \code{colGeometry}, \code{annotGeometry}, and \code{reducedDim}
#'   counterparts), and \code{spatialGraphs} are dropped because those results
#'   no longer apply after aggregation.
#' @note For developers: When debugging this function after calling
#'   \code{devtools::load_all(".")}, you may get an error that comes from S3
#'   dispatch of \code{aggregate.Vector} from the \code{S4Vectors} package. When
#'   that happens, either restart the R session, or run
#'   \code{setGeneric("aggregate", function(x, ...)
#'   standardGeneric("aggregate"))} in the console to make an S4 generic as done
#'   in the \code{terra} package to prioritize S4 dispatch.
#' @export
#' @importFrom stats aggregate
#' @concept Geometric operations
#' @examples
#' # example code
#' 
setMethod("aggregate", "SpatialFeatureExperiment", .aggregate_SFE)

# Function to make grid for multiple samples
.make_grid_samples <- function(x, sample_id, cellsize, square,
                               flat_topped) {
    bboxes <- .bbox2sf(bbox(x, sample_id = sample_id), sample_id)
    out <- lapply(st_geometry(bboxes), st_make_grid, cellsize = cellsize,
                  square = square, flat_topped = flat_topped)
    names(out) <- sample_id
    out
}
