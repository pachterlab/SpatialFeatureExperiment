#' Coordinate transforms
#'
#' `r lifecycle::badge("experimental")` These constructor functions specify
#' affine transforms for interactive use in the R session. This is used to
#' specify affine transforms in \code{\link{BioFormatsImage}} and to record
#' affine transforms performed to the SFE object as a whole.
#'
#' @name AffineTransform
#' @inheritParams rotateImg
#' @inheritParams scaleImg
#' @inheritParams translateImg
#' @inheritParams affineImg
#' @importFrom methods setOldClass
#' @param direction Character, should be either "vertical" to flip the y axis or
#'   "horizontal" to flip the x axis.
#' @concept Specification of affine transformations
#' @exportClass AffineTransform
setOldClass("AffineTransform")

#' @rdname AffineTransform
#' @export
ct_affine <- function(M, v) {
    structure(list(M = M, v = v), class = "AffineTransform")
}

#' @rdname AffineTransform
#' @export
ct_transpose <- function() {
    structure(list(name = "transpose"), class = "AffineTransform")
}

#' @rdname AffineTransform
#' @export
ct_mirror <- function(direction = c("vertical", "horizontal")) {
    direction <- match.arg(direction)
    structure(list(name = "mirror", direction = direction), class = "AffineTransform")
}

#' @rdname AffineTransform
#' @export
ct_rotate <- function(degrees) {
    structure(list(name = "rotate", degrees = degrees), class = "AffineTransform")
}

#' @rdname AffineTransform
#' @export
ct_scale <- function(factor) {
    structure(list(name = "scale", factor = factor), class = "AffineTransform")
}

#' @rdname AffineTransform
#' @export
ct_identity <- function() {
    structure(list(), class = "AffineTransform")
}

#' Combine affine transform specifications
#'
#' This function combines two affine transforms as specified in
#' \code{AffineTransform} constructors such as \code{\link{ct_affine}} into one
#' affine transform, specified by a new linear transformation matrix and
#' translation vector. The \code{bbox} argument is used to ensure that the
#' center of the bounding box specified here is preserved after transforms,
#' except when the matrix and vector are directly specified.
#'
#' @param x An \code{AffineTransform} object
#' @param y Another \code{AffineTransform} object
#' @param bbox Bounding box whose center is preserved (except for
#'   \code{ct_affine}), a numeric vector with names xmin, xmax, ymin, and ymax
#'   in any order. If NULL, then the center is assumed to be (0,0).
#' @return An \code{AffineTransform} object with the matrix and vector
#'   specifying the new combined transform.
#' @export
combineTransforms <- function(x, y, bbox) {
    if (!length(x)) return(y)
    if (!length(y)) return(x)
    old <- .get_mv(old, bbox)
    new <- .get_mv(new, bbox)
    out <- ct_affine(M = new$M %*% old$M, v = new$M %*% old$v + new$v)
    if (isTRUE(all.equal(out$M, diag(nrow = 2))) &&
        isTRUE(all.equal(as.vector(out$v), c(0,0)))) {
        out <- ct_identity()
    }
    out
}
