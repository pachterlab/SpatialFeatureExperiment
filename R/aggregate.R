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
#' @export
aggregateTx <- function(file, df = NULL, by = NULL, sample_id = "sample01",
                        spatialCoordsNames = c("X", "Y", "Z"),
                        gene_col = "gene",
                        phred_col = "qv", min_phred = 20, flip_geometry = FALSE,
                        cellsize = NULL, square = TRUE, flat_topped = FALSE,
                        new_geometry_name = "bins", unit = "micron") {
    # This is only for one file, one sample
    if (!is.null(df)) file <- df
    mols <- .check_tx_file(file, spatialCoordsNames, gene_col, phred_col,
                           min_phred, flip_geometry, BPPARAM)
    mols <- df2sf(mols, spatialCoordsNames = spatialCoordsNames,
                  geometryType = "POINT")
    if (is.null(by))
        by <- st_make_grid(mols, cellsize = cellsize, square = square,
                           flat_topped = flat_topped)
    else if (is(by, "sf")) by <- st_geometry(by)
    grid_sf <- st_sf(grid_id = seq_along(by), geometry = by)
    mols <- st_join(mols, grid_sf) # Took 5.87 minutes for 7171453 spots and 8555 bins
    mols <- st_drop_geometry(mols) |> as.data.table()
    mols <- mols[, .N, by = .(gene, grid_id)]
    mols$gene <- factor(mols$gene) # The levels are alphabetically arranged
    mols$gene_index <- as.integer(mols$gene)
    new_mat <- sparseMatrix(i = mols$gene_index, j = mols$grid_id, x = mols$N)
    rownames(new_mat) <- levels(mols$gene)
    colnames(new_mat) <- seq_len(ncol(new_mat))
    new_mat <- new_mat[,colSums(new_mat) > 0] # Remove empty grid cells
    cgs <- list(bins = grid_sf[grid_sf$grid_id %in% mols$grid_id, "geometry"])
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
                            new_geometry_name = "bins") {
    if (!is.null(df)) data_dir <- NULL
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
        if (fun_name == "sum")
            mat_agg <- sparseMatrix(i = unlist(inds), j = seq_along(inds), x = 1,
                                    dims = c(ncol(mat), length(inds)))
        else if (fun_name == "mean") {
            ll <- lengths(inds)
            mat_agg <- sparseMatrix(i = unlist(inds), j = seq_along(inds),
                                    x = rep(1/ll, times = ll),
                                    dims = c(ncol(mat), length(inds)))
        }
        out_agg <- mat %*% mat_agg
    } else if (fun_name %in% getNamespaceExports("sparseMatrixStats")) {
        # For functions supported by sparseMatrixStats
        if (!grepl("^row", fun_name))
            stop("Only row* functions from sparseMatrixStats can be used for FUN.")
        cnts_agg <- bplapply(inds, function(i) {
            m <- mat[,i,drop = FALSE]
            FUN(m)
        }, BPPARAM = BPPARAM)
        out_agg <- matrix(unlist(out_agg), ncol = length(inds))
    } else stop("Function ", fun_name, " not supported for aggregating SFE.")
    out_agg
}

# For one sample
.aggregate_sample_cell <- function(x, by, FUN, fun_name, colGeometryName, join,
                                   new_geometry_name, BPPARAM) {
    cg <- colGeometry(x, type = colGeometryName)
    inds <- join(by, cg) # somewhat faster than join(cg, by) when by has fewer geometries than cg
    not_empty <- which(length(inds))
    by <- by[not_empty]
    inds <- inds[not_empty]
    cnts <- counts(x)
    # Aggregate the gene count matrix, only do the counts assay
    # Express this as a matrix multiplication for sum. Use sweep to deal with mean
    cnts_agg <- .aggregate_num(cnts, inds, FUN, fun_name, BPPARAM)
    # Aggregate numeric columns of colData; logical are converted to numeric
    df <- colData(x)
    which_logical <- which(vapply(df, is.logical, FUN.VALUE = logical(1)))
    for (i in which_logical) {
        df[,i] <- as.integer(df[,i])
    }
    numeric_cols <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    if (any(numeric_cols)) {
        cd <- t(as.matrix(df[,numeric_cols, drop = FALSE]))
        cd_agg <- .aggregate_num(cd, inds, FUN, fun_name, BPPARAM) |> t()
        cd_agg <- as(cd_agg, "DFrame")
    }
    # Deal with anything that is neither numerical nor logical
    # Primarily character and factor. What about list columns? I don't know.
    if (!all(numeric_cols)) {
        df_bin <- data.frame(bin = rep(seq_along(inds), times = lengths(inds)),
                             index = unlist(inds))
        df_bin <- df_bin[order(df_bin$index),]
        names_not_num <- names(df)[!numeric_cols]
        for (n in names_not_num)
            df_bin[[n]] <- split(df[[,i]][n], list(bin = df_bin$bin))
        if (any(numeric_cols)) cd_agg <- cbind(df_bin, cd_agg)
    }
    cgs <- list(bins = by)
    names(cgs) <- new_geometry_name
    SpatialFeatureExperiment(assays = list(counts = cngs_agg), colData = cd_agg,
                             sample_id = sampleIDs(x),
                             colGeometries = cgs)
}

# Might turn this into an exported function
.aggregate_sample_tx <- function(x, by, rowGeometryName, new_geometry_name) {
    rg <- rowGeometry(x, rowGeometryName)
    if (!is.null(st_z_range(rg)))
        by <- st_zm(by, drop = FALSE, what = "Z")
    grid_sf <- st_sf(grid_id = seq_along(by), geometry = by)

    # Probably faster than directly calling st_intersection, since I don't need
    # the actual geometries of the intersections, maybe not
    tx_coords <- st_coordinates(rg) # Might write another function similar to formatTxTech to skip this
    tx_point <- df2sf(tx_coords, spatialCoordsNames = c("X", "Y", "Z"))
    tx_ind <- as.data.table(st_drop_geometry(tx_point))
    tx_info <- txSpots(sfe) |> st_drop_geometry()
    tx_info$L1 <- seq_along(tx_info$gene) # it has to be "gene" if it's from formatTxSpots
    tx_point <- merge(tx_point, tx_info, by = "L1") # takes a while
    tx_point <- st_as_sf(tx_point) |> st_join(grid_sf) # takes a few minutes
    tx_counts <- tx_point |> st_drop_geometry() |> as.data.table()
    tx_counts <- tx_counts[, .N, by = .(gene, L1, grid_id)]

    new_mat <- sparseMatrix(i = tx_counts$L1, j = tx_counts$grid_id, x = tx_counts$N)
    cgs <- list(bins = by)
    names(cgs) <- new_geometry_name
    SpatialFeatureExperiment(assays = list(counts = new_mat),
                             colGeometries = cgs)
}

.aggregate_sample <- function(x, by = NULL, FUN = sum,
                              colGeometryName = 1L, rowGeometryName = NULL,
                              join = st_intersects, new_geometry_name = "bins",
                              BPPARAM = SerialParam()) {
    # Here x is an SFE object with one sample
    # by is sfc
    # Can't do S4 method with signature for `by` because the argument `by` isn't
    # in the generic and I don't want to mess with the `aggregate` function in
    # other packages
    if (is.null(rowGeometryName)) {
        .aggregate_sample_tx(x, by, rowGeometryName, new_geometry_name)
    } else {
        .aggregate_sample_cell(x, by, FUN, fun_name, colGeometryName, join,
                                new_geometry_name, BPPARAM)
    }
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
#'   the gene count matrix. This can be \code{sum}, \code{mean}, or any
#'   \code{row*} function in the \code{sparseMatrixStats} package such as
#'   \code{\link{rowMedians}}. Aggregation is not done when aggregating by
#'   transcript spots in \code{rowGeometry}.
#' @param sample_id Which samples to aggregate, defaults to "all".
#' @param colGeometryName Which \code{colGeometry} to spatially aggregate the
#'   data, by default the first one.
#' @param rowGeometryName Which \code{rowGeometry} to spatially aggregate
#' @param new_geometry_name Name to give to the new \code{colGeometry} in the
#' output. Defaults to "bins".
#' @param BPPARAM A \code{\link{BiocParallelParam}} object specifying parallel
#'   computing when aggregating different genes. Defaults to
#'   \code{SerialParam()}.
#' @return An SFE object with \code{colGeometry} the same as the geometry
#'   specified in \code{by} or same as the grid specified in \code{cellsize}.
#'   \code{rowGeometries} and \code{rowData} remain the same as in the input
#'   \code{x}. \code{reducedDims}, \code{localResults}, \code{colFeatureData}
#'   (and its \code{colGeometry}, \code{annotGeometry}, and \code{reducedDim}
#'   counterparts), and \code{spatialGraphs} are dropped because those results
#'   no longer apply after aggregation.
#' @export
#' @concept Geometric operations
#' @examples
#' # example code
#'
setMethod("aggregate", "SpatialFeatureExperiment",
          function(x, by = NULL, FUN = sum, sample_id = "all",
                   colGeometryName = 1L, rowGeometryName = NULL,
                   cellsize = NULL, square = TRUE, flat_topped = FALSE,
                   new_geometry_name = "bins",
                   BPPARAM = SerialParam()) {
              sample_id <- .check_sample_id(x, sample_id, one = FALSE)
              if (is.null(by) && is.null(cellsize)) {
                  stop("Either `by` or `cellsize` must be specified.")
              }
              # Make grid for multiple samples if `by` is not specified
              if (is.null(by)) {
                  by <- .make_grid_samples(x, sample_id, colGeometryName,
                                           cellsize, square, flat_topped)
              }
              if (!is(by, "sfc") && !is(by, "sf"))
                  stop("`by` must be either sf or sfc.")
              if (length(sample_id) > 1L && (is(by, "sfc") || !"sample_id" %in% names(by))) {
                  stop("`by` must be a sf data frame with a column `sample_id`")
              }
              # Need to check new_geometry_name
              fun_name <- substitute(FUN)
              sfes <- splitSamples(x) # Output list should have sample IDs as names
              if (length(sample_id) > 1L) {
                  by <- split(by$geometry, list(sample_id = by$sample_id))
              }
              sfes <- lapply(sample_id, function(s) {
                  .aggregate_sample(sfes[[s]], by = by[[s]], FUN = FUN,
                                    colGeometryName = colGeometryName,
                                    rowGeometryName = rowGeometryName)
              })
              out <- do.call(cbind, sfes)
              # Add the original rowGeometries back
              rowGeometries(out) <- rowGeometries(x, sample_id = sample_id)
              # Keep imgData
              id_orig <- imgData(x)
              imgData(out) <- id_orig[id_orig$sample_id %in% sample_id,]
              out
          })

# Function to make grid for multiple samples
.make_grid_samples <- function(x, sample_id, colGeometryName, cellsize, square,
                               flat_topped) {
    cg <- colGeometry(x, sample_id = "all", type = colGeometryName)
    cg$sample_id <- x$sample_id
    cg <- cg[cg$sample_id %in% sample_id,]
    cgs <- split(cg, f = cg$sample_id)
    lapply(cgs, st_make_grid, cellsize = cellsize, square = square, flat_topped = flat_topped)
}
