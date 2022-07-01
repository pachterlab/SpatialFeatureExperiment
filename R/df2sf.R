.df2sf_check <- function(df, spatialCoordsNames, geometryType) {
  # For anything other than points
  if (!"ID" %in% names(df)) {
    stop("Column 'ID' must be present when specifying ", geometryType, "s.")
  }
  if (grepl("MULTI", geometryType) && !"group" %in% names(df)) {
    stop("Column 'group' must be present when specifying ", geometryType, "s.")
  }
  if (any(!names(df) %in% c("ID", "subID", "sample_id", "group", spatialCoordsNames))) {
    message("Geometry attributes are ignored.")
  }
  n_vertices <- table(df$ID)
  ids <- names(n_vertices)
  min_vertices <- switch (geometryType,
                          LINESTRING = 2L,
                          POLYGON = 3L,
                          MULTIPOINT = 1L,
                          MULTILINESTRING = 2L,
                          MULTIPOLYGON = 3L
  )
  ids_rm <- names(n_vertices[n_vertices < min_vertices])
  if (length(ids_rm)) {
    warning("Removed ", length(ids_rm), " items that have fewer than the minimum of ",
            min_vertices, " vertices for geometry type ", geometryType)
  }
  df <- df[!df$ID %in% ids_rm,]
  if (!nrow(df)) {
    stop("All geometries have fewer than ", min_vertices, " vertices. ",
         "Cannot construct ", geometryType)
  }
  df
}

.df2sf_point <- function(df, spatialCoordsNames, spotDiameter, multi) {
  # Case 1: centroids, use POINT
  if (!is.na(spotDiameter)) {
    if (spotDiameter <= 0) {
      stop("spotDiameter must be a positive number.")
    }
  }
  if (multi) {
    df <- .df2sf_check(df, spatialCoordsNames, "MULTIPOINT")
    df_split <- split(df, df$group)
    geometry_use <- lapply(df_split, function(x) {
      st_multipoint(as.matrix(x[,spatialCoordsNames]))
    })
    geometry_use <- st_sfc(geometry_use)
    out <- .df_split_sample_id(names(df), df_split, geometry_use)
  } else {
    df$geometry <- lapply(seq_len(nrow(df)), function(i) {
      st_point(c(df[[spatialCoordsNames[1]]][i], df[[spatialCoordsNames[2]]][i]))
    })
    df$geometry <- st_sfc(df$geometry)
    # Remove the original coordinate columns
    df[,spatialCoordsNames] <- NULL
    out <- st_sf(df, sf_column_name = "geometry", row.names = rownames(df))
  }
  if (!is.na(spotDiameter)) {
    out$geometry <- st_buffer(out$geometry, spotDiameter/2)
  }
  out
}

.df2poly_mat <- function(x, spatialCoordsNames) {
  m <- as.matrix(x[,spatialCoordsNames])
  # Close the polygon
  rbind(m, m[1,])
}

.df_split_sample_id <- function(nms, df_split, geometry_use) {
  if ("sample_id" %in% nms) {
    sample_ids <- vapply(df_split, function(d) unique(d$sample_id),
                         FUN.VALUE = character(1))
    out <- st_sf(ID = names(df_split), sample_id = sample_ids,
                 geometry = geometry_use,  crs = NA,
                 row.names = names(df_split))
  } else {
    out <- st_sf(ID = names(df_split), geometry = geometry_use, crs = NA,
                 row.names = names(df_split))
  }
  return(out)
}
.df2sf_polygon <- function(df, spatialCoordsNames, multi) {
  df <- unique(df)
  gt <- if (multi) "MULTIPOLYGON" else "POLYGON"
  df <- .df2sf_check(df, spatialCoordsNames, gt)
  if (multi) {
    df_split <- split(df, df$group)
    geometry_use <- lapply(df_split, function(x) {
      ms1 <- split(x, x$ID)
      if ("subID" %in% names(df)) {
        m <- lapply(ms1, function(y) {
          ms2 <- split(y, y$subID)
          lapply(ms2, .df2poly_mat, spatialCoordsNames = spatialCoordsNames)
        })
      } else {
        m <- lapply(ms1, function(y) list(.df2poly_mat(y, spatialCoordsNames)))
      }
      st_multipolygon(m)
    })
  } else {
    df_split <- split(df, df$ID)
    if ("subID" %in% names(df)) {
      # Might be holes
      geometry_use <- lapply(df_split, function(y) {
        ms <- split(y, y$subID)
        out <- lapply(ms, .df2poly_mat, spatialCoordsNames = spatialCoordsNames)
        st_polygon(out)
      })
    } else { # Definitely no holes
      geometry_use <- lapply(df_split, function(y) {
        st_polygon(list(.df2poly_mat(y, spatialCoordsNames)))
      })
    }
  }
  geometry_use <- st_sfc(geometry_use)
  .df_split_sample_id(names(df), df_split, geometry_use)
}

.df2sf_linestring <- function(df, spatialCoordsNames, multi) {
  df <- unique(df)
  gt <- if (multi) "MULTILINESTRING" else "LINESTRING"
  df <- .df2sf_check(df, spatialCoordsNames, gt)
  if (multi) {
    df_split <- split(df, df$group)
    geometry_use <- lapply(df_split, function(x) {
      ms <- split(x, x$ID)
      m <- lapply(ms, function(y) as.matrix(y[,spatialCoordsNames]))
      st_multilinestring(m)
    })
  } else {
    df_split <- split(df, df$ID)
    geometry_use <- lapply(df_split, function(x) {
      st_linestring(as.matrix(x[,spatialCoordsNames]))
    })
  }
  geometry_use <- st_sfc(geometry_use)
  .df_split_sample_id(names(df), df_split, geometry_use)
}

.is_de_facto_point <- function(df) {
  (!"ID" %in% names(df) | !anyDuplicated(df$ID)) & !"group" %in% names(df)
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
#' @param df An ordinary data frame, i.e. not \code{sf}. Or a matrix that can be
#' converted to a data frame.
#' @param spatialCoordsNames Column names in \code{df} that specify spatial
#' coordinates.
#' @param geometryType Type of geometry to convert the ordinary data frame to.
#' If the geometry in \code{df} is de facto points, then this argument will be
#' ignored and the returned \code{sf} will have geometry type POINT. For any
#' geometry type where one geometry is specified by multiple coordinates, the
#' data frame \code{df} must have a column "ID" specifying which coordinate
#' belongs to which geometry. For MULTI* geometries, there must be a "group"
#' column specifying which coordinates for which MULTI geometry.
#' @return An \code{sf} object.
#' @export
#' @examples
#' # Points, use spotDiameter to convert to circle polygons
#' # This is done to Visium spots
#' pts_df <- readRDS(system.file("testdata/pts_df.rds",
#'                   package = "SpatialFeatureExperiment"))
#' sf_use <- df2sf(pts_df, geometryType = "POINT", spotDiameter = 0.1)
#' # Linestring
#' ls_df <- readRDS(system.file("testdata/ls_df.rds",
#'                  package = "SpatialFeatureExperiment"))
#' sf_use <- df2sf(ls_df, geometryType = "LINESTRING")
#' # Polygon
#' pol_df <- readRDS(system.file("testdata/pol_df.rds",
#'                   package = "SpatialFeatureExperiment"))
#' sf_use <- df2sf(pol_df, geometryType = "POLYGON",
#'                 spatialCoordsNames = c("V1", "V2"))
#' # Multipolygon
#' mpol_df <- readRDS(system.file("testdata/mpol_df.rds",
#'                    package = "SpatialFeatureExperiment"))
#' sf_use <- df2sf(mpol_df, geometryType = "MULTIPOLYGON",
#'                 spatialCoordsNames = c("V1", "V2"))
#' # Multiple sample_ids present
#' multipts_df <- readRDS(system.file("testdata/multipts_df.rds",
#'                        package = "SpatialFeatureExperiment"))
#' sf_use <- df2sf(multipts_df, geometryType = "MULTIPOINT")
df2sf <- function(df, spatialCoordsNames = c("x", "y"), spotDiameter = NA,
                  geometryType = c("POINT", "LINESTRING", "POLYGON",
                                   "MULTIPOINT", "MULTILINESTRING",
                                   "MULTIPOLYGON")) {
  if (is.matrix(df)) df <- as.data.frame(df)
  if (any(!spatialCoordsNames %in% names(df))) {
    cols_absent <- setdiff(spatialCoordsNames, names(df))
    if (length(cols_absent) > 1L) {
      stop("Columns ", paste(cols_absent, collapse = ", "), " are absent.")
    } else {
      stop("Column ", cols_absent, " is absent.")
    }
  }
  if (.is_de_facto_point(df)) geometryType <- "POINT"
  geometryType <- match.arg(geometryType)
  switch (geometryType,
          POINT = .df2sf_point(df, spatialCoordsNames, spotDiameter, multi = FALSE),
          MULTIPOINT = .df2sf_point(df, spatialCoordsNames, spotDiameter, multi = TRUE),
          LINESTRING = .df2sf_linestring(df, spatialCoordsNames, multi = FALSE),
          MULTILINESTRING = .df2sf_linestring(df, spatialCoordsNames, multi = TRUE),
          POLYGON = .df2sf_polygon(df, spatialCoordsNames, multi = FALSE),
          MULTIPOLYGON = .df2sf_polygon(df, spatialCoordsNames, multi = TRUE)
  )
}

# Call in SFE constructor and *Geometries replacement methods
.df2sf_in_list <- function(x, spatialCoordsNames = c("x", "y"),
                           spotDiameter = NA, geometryType = "POLYGON") {
  if (!is.null(x) && !is(x, "sf") && !is.data.frame(x)) {
    stop("Each element of the list for *Geometry must be an ",
         "sf object or a data frame.")
  }
  if (is(x, "sf") || is.null(x)) {
    return(x)
  } else if (is.data.frame(x)) {
    return(df2sf(x, spatialCoordsNames, spotDiameter, geometryType))
  }
}

.df2sf_list <- function(x, spatialCoordsNames = c("x", "y"),
                        spotDiameter = NA, geometryType = "POLYGON") {
  x_is_sf <- vapply(x, function(t) is(t, "sf"), FUN.VALUE = logical(1))
  if (all(x_is_sf)) {
    return(x)
  }
  if (length(geometryType) == 1L) {
    geometryType <- rep(geometryType, length(x))
  } else if (length(geometryType) != length(x)) {
    stop("geometryTypes must be either length 1 or the same length ",
         "as the input list.")
  }
  mapply(.df2sf_in_list, x = x, geometryType = geometryType,
         MoreArgs = list(spatialCoordsNames = spatialCoordsNames,
                         spotDiameter = spotDiameter),
         SIMPLIFY = FALSE)
}
