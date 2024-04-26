.df2sf_check <- function(df, spatialCoordsNames, geometryType,
                         group_col = "group", id_col = "ID", subid_col = "subID") {
    # For anything other than points
    if (geometryType == "MULTIPOINT") {
        id_col <- group_col # should be the same
    } else if (!id_col %in% names(df)) {
        stop("Column ", id_col, " for individual geometries is absent.")
    }
    if (grepl("MULTI", geometryType) && !group_col %in% names(df)) {
        stop("Column", group_col, " to identify MULTI geometries is abesent.")
    }
    df <- unique(df)
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
    # Only keep other attributes with one value per geometry
    cols_keep <- c(group_col, id_col, subid_col,
                   "sample_id", spatialCoordsNames)
    cols_check <- setdiff(names(df), c(group_col, id_col, subid_col,
                                       "sample_id", spatialCoordsNames))
    col_geo <- if (group_col %in% names(df)) group_col else id_col
    n_geos <- length(unique(df[[col_geo]]))
    if_keep <- vapply(cols_check,
                      function(x) {
                          if (is.data.table(df))
                              df_check <- df[,c(x, col_geo), with=FALSE]
                          else df_check <- df[,c(x, col_geo)]
                          nrow(unique(df_check)) == n_geos
                      }, FUN.VALUE = logical(1))
    cols_use <- c(cols_check[if_keep], intersect(names(df), cols_keep))
    # To work around data.table's nicer column subsetting with symbols
    # and to remain compatible with base data frames
    if (!is.data.table(df)) ..cols_use <- cols_use
    df[,..cols_use]
}

.use_z <- function(spatialCoordsNames) {
    if (length(spatialCoordsNames) > 2L) spatialCoordsNames[3] else NULL
}

.use_subid <- function(df, subid_col) {
    if (subid_col %in% names(df)) subid_col else NULL
}

.df2sf_point <- function(df, spatialCoordsNames, spotDiameter, multi,
                         group_col = "group", check = TRUE, ...) {
    # Case 1: centroids, use POINT
    if (!is.na(spotDiameter)) {
        if (spotDiameter <= 0) {
            stop("spotDiameter must be a positive number.")
        }
    }
    if (multi) {
        if (check) df <- .df2sf_check(df, spatialCoordsNames, "MULTIPOINT",
                                      group_col = group_col)
        if (!is.data.table(df)) ..group_col <- group_col
        df <- df[order(df[,..group_col]),]
        out <- sf_multipoint(df, x = spatialCoordsNames[1],
                             y = spatialCoordsNames[2],
                             z = .use_z(spatialCoordsNames),
                             multipoint_id = group_col,
                             keep = TRUE)
        rownames(out) <- out[[group_col]]
    } else {
        rns <- rownames(df)
        out <- sf::st_as_sf(df, coords = spatialCoordsNames, crs = NA,
                            row.names = rns) # in this case faster than sf_point
    }
    if (!is.na(spotDiameter)) {
      st_geometry(out) <- st_buffer(st_geometry(out), spotDiameter / 2, ...)
    }
    out
}

.df2sf_polygon <- function(df, spatialCoordsNames, multi,
                           group_col = "group", id_col, subid_col,
                           check = TRUE) {
    gt <- if (multi) "MULTIPOLYGON" else "POLYGON"
    if (check) df <- .df2sf_check(df, spatialCoordsNames, gt,
                                  group_col, id_col, subid_col)
    if (multi) {
        out <- sf_multipolygon(df, x = spatialCoordsNames[1],
                               y = spatialCoordsNames[2],
                               z = .use_z(spatialCoordsNames),
                               polygon_id = id_col,
                               multipolygon_id = group_col,
                               linestring_id = .use_subid(df, subid_col),
                               keep = TRUE)
        rownames(out) <- out[[group_col]]
    } else {
        out <- sf_polygon(df, x = spatialCoordsNames[1],
                          y = spatialCoordsNames[2],
                          z = .use_z(spatialCoordsNames),
                          polygon_id = id_col,
                          linestring_id = .use_subid(df, subid_col),
                          keep = TRUE)
        rownames(out) <- out[[id_col]]
    }
    return(out)
}

.df2sf_linestring <- function(df, spatialCoordsNames, multi,
                              group_col = "group", id_col, check = TRUE) {
    gt <- if (multi) "MULTILINESTRING" else "LINESTRING"
    if (check) df <- .df2sf_check(df, spatialCoordsNames, gt,
                                  group_col, id_col)
    if (multi) {
        out <- sf_multilinestring(df, x = spatialCoordsNames[1],
                                  y = spatialCoordsNames[2],
                                  z = .use_z(spatialCoordsNames),
                                  multilinestring_id = group_col,
                                  linestring_id = id_col,
                                  keep = TRUE)
        rownames(out) <- out[[group_col]]
    } else {
        out <- sf_linestring(df, x = spatialCoordsNames[1],
                             y = spatialCoordsNames[2],
                             z = .use_z(spatialCoordsNames),
                             linestring_id = id_col,
                             keep = TRUE)
        rownames(out) <- out[[id_col]]
    }
    return(out)
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
#'   converted to a data frame.
#' @param spatialCoordsNames Column names in \code{df} that specify spatial
#'   coordinates.
#' @param geometryType Type of geometry to convert the ordinary data frame to.
#'   If the geometry in \code{df} is de facto points, then this argument will be
#'   ignored and the returned \code{sf} will have geometry type POINT.
#' @param group_col Column to indicate which coordinates for which MULTI
#'   geometry, such as to identify which MULTIPOLYGON or MULTIPOINT.
#' @param id_col Column to indicate coordinates for which geometry, within a
#'   MULTI geometry if applicable, such as to identify which POLYGON or which
#'   polygon within a MULTIPOLYGON.
#' @param check Logical, whether to check the input data frame for issues
#'   related to constructing the geometry of interese such as number of vertices
#'   per geometry. If \code{FALSE}, it will save a bit of time, which is useful
#'   when the input is already known to be good.
#' @param subid_col Column to indicate coordinates for holes in polygons.
#' @param ... Other arguments passed to `sf::st_buffer`, mainly to make polygon shapes,
#'  eg Visium spot `endCapStyle = "ROUND"` and VisiumHD bin `endCapStyle = "SQUARE"`
#' @return An \code{sf} object.
#' @export
#' @concept Utilities
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom sfheaders sf_multipoint sf_polygon sf_multipolygon sf_linestring
#'   sf_multilinestring
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
                  ),
                  group_col = "group",
                  id_col = "ID",
                  subid_col = "subID", check = TRUE,
                  BPPARAM = deprecated(), ...) {
    if (is_present(BPPARAM)) {
        deprecate_warn("1.6.0", "df2sf(BPPARAM = )",
                       details = "The sfheaders package is now used instead for much better performance")
    }
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
    out <- switch(geometryType,
        POINT = .df2sf_point(df, spatialCoordsNames, spotDiameter, multi = FALSE, ...),
        MULTIPOINT = .df2sf_point(df, spatialCoordsNames, spotDiameter, multi = TRUE,
                                  group_col = group_col, check = check, ...),
        LINESTRING = .df2sf_linestring(df, spatialCoordsNames, multi = FALSE,
                                       id_col = id_col, check = check),
        MULTILINESTRING = .df2sf_linestring(df, spatialCoordsNames, multi = TRUE,
                                            group_col = group_col, id_col = id_col,
                                            check = check),
        POLYGON = .df2sf_polygon(df, spatialCoordsNames, multi = FALSE,
                                 id_col = id_col, subid_col = subid_col, check = check),
        MULTIPOLYGON = .df2sf_polygon(df, spatialCoordsNames, multi = TRUE,
                                      group_col = group_col, id_col = id_col,
                                      subid_col = subid_col, check = check)
    )
    out
}

# Call in SFE constructor and *Geometries replacement methods
.df2sf_in_list <- function(x, spatialCoordsNames = c("x", "y"),
                           spotDiameter = NA, geometryType = "POLYGON",
                           group_col = "group", id_col = "ID", subid_col = "subID",
                           check = TRUE) {
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
                     group_col, id_col, subid_col, check = check))
    }
}

.df2sf_list <- function(x, spatialCoordsNames = c("x", "y"),
                        spotDiameter = NA, geometryType = "POLYGON",
                        group_col = "group", id_col = "ID", subid_col = "subID",
                        check = TRUE) {
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
            group_col = group_col, id_col = id_col, subid_col = subid_col,
            check = check
        ),
        SIMPLIFY = FALSE
    )
}
