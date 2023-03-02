#' Update a SpatialFeatureExperiment object
#'
#' Update a \linkS4class{SpatialFeatureExperiment} to the latest version of
#' object structure. This is usually called by internal functions.
#'
#' Version 1.1.4 adds package version to the SFE object. We are considering an
#' overhaul of the \code{spatialGraphs} slot in a future version using the
#' \code{sfdep} package to decouple the adjacency graph from the edge weights.
#'
#' @param object An old \linkS4class{SpatialFeatureExperiment} object.
#' @param ... Additional arguments that are ignored.
#' @param verbose Logical scalar indicating whether a message should be emitted
#'   as the object is updated.
#' @return An updated version of \code{object}.
#' @seealso \code{\link{objectVersion}}, which is used to determine if the
#'   object is up-to-date.
#' @importFrom BiocGenerics updateObject
#' @export
#' @name updateObject
#' @aliases updateObject,SpatialFeatureExperiment-method
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData("small")
#' # First version of SFE object doesn't log SFE package version, so should be NULL
#' SFEVersion(sfe)
#' sfe <- updateObject(sfe)
#' # See current version
#' SFEVersion(sfe)
setMethod("updateObject", "SpatialFeatureExperiment",
          function(object, ..., verbose = FALSE) {
              old_version <- int_metadata(object)$SFE_version
              triggered <- FALSE
              curr_version <- packageVersion("SpatialFeatureExperiment")
              if (is.null(old_version)) {
                  # Meaning prior to 1.1.4, i.e. the first version
                  old_version <- package_version("1.0.0")
                  #int_metadata(object)$SFE_version <- curr_version
                  triggered <- TRUE
              }
              if (verbose && triggered) {
                  message("[updateObject] ", class(object)[1], " object uses ",
                          "internal representation\n",
                          "[updateObject] from SpatialFeatureExperiment ",
                          old_version, ". ", "Updating it to version ",
                          curr_version, "\n", appendLF = FALSE)
              }
              int_metadata(object)$SFE_version <- curr_version
              callNextMethod()
              #object
          })

#' @rdname updateObject
#' @export
SFEVersion <- function(object) int_metadata(object)$SFE_version
