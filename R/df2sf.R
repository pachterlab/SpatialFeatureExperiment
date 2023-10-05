.df2sf_check <- function(df, spatialCoordsNames, geometryType,
                         group_col = "group", id_col = "ID", subid_col = "subID") {
    # For anything other than points
    if (!id_col %in% names(df) && geometryType != "MULTIPOINT") {
        stop("Column ", id_col, " for individual geometries is absent.")
    }
    if (grepl("MULTI", geometryType) && !group_col %in% names(df)) {
        stop("Column", group_col, " to identify MULTI geometries is abesent.")
    }
    n_vertices <- table(df[[id_col]])
    ids <- names(n_vertices)
    min_vertices <- switch(geometryType,
        LINESTRING = 2L,
        POLYGON = 3L,
        MULTIPOINT = 1L,
        MULTILINESTRING = 2L,
        MULTIPOLYGON = 3L
    )
    ids_rm <- names(n_vertices[n_vertices < min_vertices])
    if (length(ids_rm)) {
        warning(
            "Removed ", length(ids_rm),
            " items that have fewer than the minimum of ",
            min_vertices, " vertices for geometry type ", geometryType
        )
    }
    df <- df[!df[[id_col]] %in% ids_rm, ]
    if (!nrow(df)) {
        stop(
            "All geometries have fewer than ", min_vertices, " vertices. ",
            "Cannot construct ", geometryType
        )
    }
    df
}

.df2sf_point <- function(df, spatialCoordsNames, spotDiameter, multi,
                         BPPARAM, group_col = "group") {
    # Case 1: centroids, use POINT
    if (!is.na(spotDiameter)) {
        if (spotDiameter <= 0) {
            stop("spotDiameter must be a positive number.")
        }
    }
    if (multi) {
        df <- .df2sf_check(df, spatialCoordsNames, "MULTIPOINT",
                           group_col = group_col)
        df_split <- split(df, df[[group_col]])
        geometry_use <- bplapply(df_split, function(x) {
            st_multipoint(as.matrix(x[, spatialCoordsNames]))
        }, BPPARAM = BPPARAM)
        geometry_use <- st_sfc(geometry_use)
        out <- .df_split_sample_id(names(df), df_split, geometry_use, spatialCoordsNames,
                                   group_col = group_col)
    } else {
        df$geometry <- bplapply(seq_len(nrow(df)), function(i) {
            st_point(unlist(df[i, spatialCoordsNames], use.names = FALSE))
        }, BPPARAM = BPPARAM)
        df$geometry <- st_sfc(df$geometry)
        # Remove the original coordinate columns
        df[, spatialCoordsNames] <- NULL
        out <- st_sf(df, sf_column_name = "geometry", row.names = rownames(df))
    }
    if (!is.na(spotDiameter)) {
        out$geometry <- st_buffer(out$geometry, spotDiameter / 2)
    }
    out
}

.df2poly_mat <- function(x, spatialCoordsNames) {
    m <- as.matrix(x[, spatialCoordsNames])
    if (!isTRUE(all.equal(m[1,], m[nrow(m)]))) {
        # Close the polygon
        rbind(m, m[1, ])
    } else m
}

.df_split_sample_id <- function(nms, df_split, geometry_use, spatialCoordsNames,
                                group_col = "group", id_col = "ID", subid_col = "subID") {
    if ("sample_id" %in% nms) {
        sample_ids <- vapply(df_split, function(d) unique(d$sample_id),
            FUN.VALUE = character(1)
        )
        out <- st_sf(
            ID = names(df_split), sample_id = sample_ids,
            geometry = geometry_use, crs = NA,
            row.names = names(df_split)
        )
    } else {
        out <- st_sf(
            ID = names(df_split), geometry = geometry_use, crs = NA,
            row.names = names(df_split)
        )
    }
    # The other attributes
    cols_use <- setdiff(names(df_split[[1]]), c(group_col, id_col, subid_col,
                                                "sample_id", spatialCoordsNames))
    if (length(cols_use)) {
        df_attrs <- lapply(df_split, function(x) {
            o <- unique(x[, cols_use, drop = FALSE])
            if (nrow(o) > 1L) {
                warning("Multiple combinations of attributes per geometry found. ",
                        "Only using the first.")
                o <- o[1,,drop = FALSE]
            }
            o
        })
        df_attrs <- do.call(rbind, df_attrs)
        out <- cbind(out, df_attrs)
    }
    return(out)
}
.df2sf_polygon <- function(df, spatialCoordsNames, multi, BPPARAM,
                           group_col = "group", id_col, subid_col) {
    df <- unique(df)
    gt <- if (multi) "MULTIPOLYGON" else "POLYGON"
    df <- .df2sf_check(df, spatialCoordsNames, gt,
                       group_col, id_col, subid_col)
    if (multi) {
        df_split <- split(df, df[[group_col]])
        geometry_use <- lapply(df_split, function(x) {
            ms1 <- split(x, x$ID)
            if ("subID" %in% names(df)) {
                m <- bplapply(ms1, function(y) {
                    ms2 <- split(y, y[[subid_col]])
                    lapply(ms2, .df2poly_mat, spatialCoordsNames = spatialCoordsNames)
                }, BPPARAM = BPPARAM)
            } else {
                m <- bplapply(ms1, function(y) list(.df2poly_mat(y, spatialCoordsNames)),
                              BPPARAM = BPPARAM)
            }
            st_multipolygon(m)
        })
    } else {
        df_split <- split(df, df[[id_col]])
        if ("subID" %in% names(df)) {
            # Might be holes
            geometry_use <- bplapply(df_split, function(y) {
                ms <- split(y, y[[subid_col]])
                out <- lapply(ms, .df2poly_mat, spatialCoordsNames = spatialCoordsNames)
                st_polygon(out)
            }, BPPARAM = BPPARAM)
        } else { # Definitely no holes
            geometry_use <- bplapply(df_split, function(y) {
                st_polygon(list(.df2poly_mat(y, spatialCoordsNames)))
            }, BPPARAM = BPPARAM)
        }
    }
    geometry_use <- st_sfc(geometry_use)
    .df_split_sample_id(names(df), df_split, geometry_use, spatialCoordsNames,
                        group_col = group_col, id_col = id_col, subid_col = subid_col)
}

.df2sf_linestring <- function(df, spatialCoordsNames, multi, BPPARAM,
                              group_col = "group", id_col) {
    df <- unique(df)
    gt <- if (multi) "MULTILINESTRING" else "LINESTRING"
    df <- .df2sf_check(df, spatialCoordsNames, gt,
                       group_col, id_col)
    if (multi) {
        df_split <- split(df, df[[group_col]])
        geometry_use <- lapply(df_split, function(x) {
            ms <- split(x, x$ID)
            m <- bplapply(ms, function(y) as.matrix(y[, spatialCoordsNames]),
                          BPPARAM = BPPARAM)
            st_multilinestring(m)
        })
    } else {
        df_split <- split(df, df[[id_col]])
        geometry_use <- bplapply(df_split, function(x) {
            st_linestring(as.matrix(x[, spatialCoordsNames]))
        }, BPPARAM = BPPARAM)
    }
    geometry_use <- st_sfc(geometry_use)
    .df_split_sample_id(names(df), df_split, geometry_use, spatialCoordsNames,
                        group_col = group_col, id_col = id_col)
}

.is_de_facto_point <- function(df, group_col, id_col) {
    (!id_col %in% names(df) || !anyDuplicated(df[[id_col]])) && !group_col %in% names(df)
}

#' From ordinary data frame to sf to construct SFE object
#'
#' While the \code{SpatialFeatureExperiment} constructor and \code{*Geometry}
#' replacement methods can convert properly formatted ordinary data frames into
#' \code{sf} objects which are used to store the geometries internally, the user
#' might want to do the conversion, check if the geometry is valid, and inspect
#' and fix any invalid geometries.
#'
#' @inheritParams SpatialFeatureExperiment
#' @inheritParams BiocParallel::bplapply
#' @param df An ordinary data frame, i.e. not \code{sf}. Or a matrix that can be
#' converted to a data frame.
#' @param spatialCoordsNames Column names in \code{df} that specify spatial
#' coordinates.
#' @param geometryType Type of geometry to convert the ordinary data frame to.
#' If the geometry in \code{df} is de facto points, then this argument will be
#' ignored and the returned \code{sf} will have geometry type POINT.
#' @param group_col Column to indicate which coordinates for which MULTI geometry,
#' such as to identify which MULTIPOLYGON or MULTIPOINT.
#' @param id_col Column to indicate coordinates for which geometry, within a
#' MULTI geometry if applicable, such as to identify which POLYGON or which
#' polygon within a MULTIPOLYGON.
#' @param subid_col Column to indicate coordinates for holes in polygons.
#' @return An \code{sf} object.
#' @export
#' @concept Utilities
#' @importFrom BiocParallel bplapply SerialParam
#' @examples
#' # Points, use spotDiameter to convert to circle polygons
#' # This is done to Visium spots
#' pts_df <- readRDS(system.file("extdata/pts_df.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' sf_use <- df2sf(pts_df, geometryType = "POINT", spotDiameter = 0.1)
#' # Linestring
#' ls_df <- readRDS(system.file("extdata/ls_df.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' sf_use <- df2sf(ls_df, geometryType = "LINESTRING")
#' # Polygon
#' pol_df <- readRDS(system.file("extdata/pol_df.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' sf_use <- df2sf(pol_df,
#'     geometryType = "POLYGON",
#'     spatialCoordsNames = c("V1", "V2")
#' )
#' # Multipolygon
#' mpol_df <- readRDS(system.file("extdata/mpol_df.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' sf_use <- df2sf(mpol_df,
#'     geometryType = "MULTIPOLYGON",
#'     spatialCoordsNames = c("V1", "V2")
#' )
#' # Multiple sample_ids present
#' multipts_df <- readRDS(system.file("extdata/multipts_df.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' sf_use <- df2sf(multipts_df, geometryType = "MULTIPOINT")
df2sf <- function(df, spatialCoordsNames = c("x", "y"), spotDiameter = NA,
                  geometryType = c(
                      "POINT", "LINESTRING", "POLYGON",
                      "MULTIPOINT", "MULTILINESTRING",
                      "MULTIPOLYGON"
                  ), BPPARAM = SerialParam(),
                  group_col = "group",
                  id_col = "ID",
                  subid_col = "subID") {
    if (is.matrix(df)) df <- as.data.frame(df)
    if (any(!spatialCoordsNames %in% names(df))) {
        cols_absent <- setdiff(spatialCoordsNames, names(df))
        if (length(cols_absent) > 1L) {
            stop("Columns ", paste(cols_absent, collapse = ", "), " are absent.")
        } else {
            stop("Column ", cols_absent, " is absent.")
        }
    }
    if (.is_de_facto_point(df, group_col, id_col)) geometryType <- "POINT"
    geometryType <- match.arg(geometryType)
    switch(geometryType,
        POINT = .df2sf_point(df, spatialCoordsNames, spotDiameter, multi = FALSE, BPPARAM),
        MULTIPOINT = .df2sf_point(df, spatialCoordsNames, spotDiameter, multi = TRUE, BPPARAM,
                                  group_col = group_col),
        LINESTRING = .df2sf_linestring(df, spatialCoordsNames, multi = FALSE, BPPARAM,
                                       id_col = id_col),
        MULTILINESTRING = .df2sf_linestring(df, spatialCoordsNames, multi = TRUE, BPPARAM,
                                            group_col = group_col, id_col = id_col),
        POLYGON = .df2sf_polygon(df, spatialCoordsNames, multi = FALSE, BPPARAM,
                                 id_col = id_col, subid_col = subid_col),
        MULTIPOLYGON = .df2sf_polygon(df, spatialCoordsNames, multi = TRUE, BPPARAM,
                                      group_col = group_col, id_col = id_col,
                                      subid_col = subid_col)
    )
}

# Call in SFE constructor and *Geometries replacement methods
.df2sf_in_list <- function(x, spatialCoordsNames = c("x", "y"),
                           spotDiameter = NA, geometryType = "POLYGON",
                           BPPARAM = SerialParam(),
                           group_col = "group", id_col = "ID", subid_col = "subID") {
    if (!is.null(x) && !is(x, "sf") && !is.data.frame(x) && !is.matrix(x)) {
        stop(
            "Each element of the list for *Geometry must be an ",
            "sf object or a data frame or a matrix."
        )
    }
    if (is(x, "sf") || is.null(x)) {
        return(x)
    } else if (is.data.frame(x) || is.matrix(x)) {
        return(df2sf(x, spatialCoordsNames, spotDiameter, geometryType,
                     BPPARAM = BPPARAM,
                     group_col, id_col, subid_col))
    }
}

.df2sf_list <- function(x, spatialCoordsNames = c("x", "y"),
                        spotDiameter = NA, geometryType = "POLYGON",
                        BPPARAM = SerialParam(),
                        group_col = "group", id_col = "ID", subid_col = "subID") {
    x_is_sf <- vapply(x, function(t) is(t, "sf"), FUN.VALUE = logical(1))
    if (all(x_is_sf)) {
        return(x)
    }
    if (length(geometryType) == 1L) {
        geometryType <- rep(geometryType, length(x))
    } else if (length(geometryType) != length(x)) {
        stop(
            "geometryTypes must be either length 1 or the same length ",
            "as the input list."
        )
    }
    mapply(.df2sf_in_list,
        x = x, geometryType = geometryType,
        MoreArgs = list(
            spatialCoordsNames = spatialCoordsNames,
            spotDiameter = spotDiameter,
            BPPARAM = BPPARAM,
            group_col, id_col, subid_col
        ),
        SIMPLIFY = FALSE
    )
}
