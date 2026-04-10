#' Rotate the tissue to reduce empty space in plots
#'
#' This function finds the minimum bounding rotated rectangle and then rotates
#' the sample so the minimum bounding rotated rectangle is horizontal or
#' vertical. This will reduce empty space in plots when the sample has a long
#' shape.
#' 
#' @inheritParams SFE-transform
#' @param x Object of interest to be rotated, can be a matrix with 2 columns,
#'   \code{sf}, \code{sfc}, or \code{SpatialFeatureExperiment}.
#' @param orientation Whether you want the rotated sample to be horizontal or
#' vertical.
#' @param ... Ignored
#' @return Object of the same class as \code{x}, rotated.
#' @importFrom sf st_multipoint st_minimum_rotated_rectangle
#' @name rotateMinRect
NULL

#' @rdname rotateMinRect
#' @export
setMethod("rotateMinRect", "matrix",
          function(x, orientation = c("horizontal", "vertical"), ...) {
              orientation <- match.arg(orientation)
              x <- x[,1:2]
              g <- st_multipoint(x) |> st_sfc()
              g <- rotateMinRect(g, orientation)
              st_coordinates(g)[,1:2]
          })

.get_theta <- function(x, orientation) {
    rect <- st_minimum_rotated_rectangle(st_union(x))
    rc <- st_coordinates(rect)
    edge_lengths <- vapply(1:2, function(i)
        (rc[i+1,1:2] - rc[i,1:2])^2 |> sum() |> sqrt(),
        FUN.VALUE = numeric(1))
    which_long <- which.max(edge_lengths)
    lc <- rc[which_long+1, 1:2] - rc[which_long, 1:2]
    if (orientation == "horizontal") {
        theta <- atan(lc[2]/lc[1])
    } else theta <- -atan(lc[1]/lc[2])
    theta
}

.rotate_sf <- function(x, orientation = c("horizontal", "vertical"), ...) {
    orientation <- match.arg(orientation)
    theta <- .get_theta(x, orientation)
    .rotate_geometry(x, bbox = st_bbox(x), degrees = theta/pi*180)
}

#' @rdname rotateMinRect
#' @export
setMethod("rotateMinRect", "sfc", .rotate_sf)

#' @rdname rotateMinRect
#' @export
setMethod("rotateMinRect", "sf", .rotate_sf)

#' @rdname rotateMinRect
#' @export
setMethod("rotateMinRect", "SpatialFeatureExperiment", 
          function(x, orientation = c("horizontal", "vertical"),
                   maxcell = 1e7) {
              orientation <- match.arg(orientation)
              sc <- spatialCoords(x)
              g <- st_multipoint(sc) |> st_sfc()
              theta <- .get_theta(g, orientation)
              SpatialFeatureExperiment::rotate(x, degrees = theta/pi*180)
          })
