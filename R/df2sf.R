.df2sf_check <- function(df, spatialCoordsNames, geometryType) {
  # For anything other than points
  if (!"ID" %in% names(df)) {
    stop("Column 'ID' must be present when specifying ", geometryType, "s.")
  }
  if (grepl("MULTI", geometryType) && !"group" %in% names(df)) {
    stop("Column 'group' must be present when specifying ", geometryType, "s.")
  }
  if (any(!names(df) %in% c("ID", spatialCoordsNames))) {
    warning("Geometry attributes are ignored.")
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
  df <- df[!df$ID %in% ids_rm,]
  if (!nrow(df)) {
    stop("All geometries have fewer than ", min_vertices, " vertices. ",
         "Cannot construct ", geometryType)
  }
}

.df2sf_point <- function(df, spatialCoordsNames, spotDiameter, multi) {
  # Case 1: centroids, use POINT
  if (!is.na(spotDiameter)) {
    if (spotDiameter <= 0) {
      stop("spotDiameter must be a positive number.")
    }
  }
  if (multi) {
    .df2sf_check(df, spatialCoordsNames, "MULTIPOINT")
    df_split <- split(df, df$group)
    geometry_use <- lapply(df_split, function(x) {
      st_multipoint(as.matrix(x[,spatialCoordsNames]))
    })
    geometry_use <- st_sfc(geometry_use)
    out <- st_sf(group = names(df_split), geometry = geometry_use)
  } else {
    df$geometry <- lapply(seq_len(nrow(df)), function(i) {
      st_point(c(df[[spatialCoordsNames[1]]], df[[spatialCoordsNames[2]]]))
    })
    df$geometry <- st_sfc(df$geometry)
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
.df2sf_polygon <- function(df, spatialCoordsNames, multi) {
  df <- unique(df)
  gt <- if (multi) "MULTIPOLYGON" else "POLYGON"
  .df2sf_check(df, spatialCoordsNames, gt)
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
  st_sf(ID = names(df_split), geometry = geometry_use, crs = NA,
        row.names = names(df_split))
}

.df2sf_linestring <- function(df, spatialCoordsNames, multi) {
  df <- unique(df)
  gt <- if (multi) "MULTILINESTRING" else "LINESTRING"
  .df2sf_check(df, spatialCoordsNames, gt)
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
  st_sf(ID = names(df_split), geometry = geometry_use, crs = NA,
        row.names = names(df_split))
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
#' @param df An ordinary data frame, i.e. not \code{sf}.
#' @param spatialCoordNames Column names in \code{df} that specify spatial
#' coordinates.
#' @param geometryType Type of geometry to convert the ordinary data frame to.
#' If the geometry in \code{df} is de facto points, then this argument will be
#' ignored and the returned \code{sf} will have geometry type POINT.
#' @return An \code{sf} object.
#' @export
df2sf <- function(df, spatialCoordsNames = c("x", "y"), spotDiameter = NA,
                  geometryType = c("POINT", "LINESTRING", "POLYGON",
                                   "MULTIPOINT", "MULTILINESTRING",
                                   "MULTIPOLYGON")) {
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
          POINT = .df2sf_point(df, spatialCoordNames, spotDiameter, multi = FALSE),
          MULTIPOINT = .df2sf_point(df, spatialCoordNames, spotDiameter, multi = TRUE),
          LINESTRING = .df2sf_linestring(df, spatialCoordNames, multi = FALSE),
          MULTILINESTRING = .df2sf_linestring(df, spatialCoordNames, multi = TRUE),
          POLYGON = .df2sf_polygon(df, spatialCoordNames, multi = FALSE),
          MULTIPOLYGON = .df2sf_polygon(df, spatialCoordNames, multi = TRUE)
  )
}

# Call in SFE constructor and *Geometries replacement methods
.df2sf_in_list <- function(x, spatialCoordsNames, spotDiameter, geometryType) {
  if (!is(x, "sf") && !is.data.frame(x)) {
    stop("Each element of the list for *Geometry must be an ",
         "sf object or a data frame.")
  }
  if (is(x, "sf")) {
    return(x)
  } else if (is.data.frame(x)) {
    return(df2sf(x, spatialCoordsNames, spotDiameter, geometryType))
  }
}

.df2sf_list <- function(x, spatialCoordsNames, spotDiameter, geometryType) {
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
  mapply(.df2sf_in_list, df = x, geometryType = geometryType,
         MoreArgs = list(spatialCoordsNames = spatialCoordsNames,
                         spotDiameter = spotDiameter))
}
