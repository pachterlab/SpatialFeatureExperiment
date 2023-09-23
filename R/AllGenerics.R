# Geometry getters and setters should be kind of similar to implementation of
# reducedDims in SCE

# Getters and setters for dimGeometry--------
# For colGeometry and rowGeometry
#' @export
setGeneric(
    "dimGeometry",
    function(x, type = 1L, MARGIN, sample_id = NULL,
             withDimnames = TRUE) {
        standardGeneric("dimGeometry")
    }
)

#' @export
setGeneric("dimGeometry<-", function(x, type = 1L, MARGIN, sample_id = NULL,
                                     withDimnames = TRUE, translate = TRUE, ...,
                                     value) {
    standardGeneric("dimGeometry<-")
})

#' @export
setGeneric(
    "dimGeometryNames",
    function(x, MARGIN) standardGeneric("dimGeometryNames")
)

#' @export
setGeneric(
    "dimGeometryNames<-",
    function(x, MARGIN, value) standardGeneric("dimGeometryNames<-")
)

#' @export
setGeneric(
    "dimGeometries",
    function(x, MARGIN, withDimnames = TRUE) {
        standardGeneric("dimGeometries")
    }
)

#' @export
setGeneric(
    "dimGeometries<-",
    function(x, MARGIN, withDimnames = TRUE, translate = TRUE, ..., value) {
        standardGeneric("dimGeometries<-")
    }
)

# Getters and setters for annotGeometry---------
#' @export
setGeneric(
    "annotGeometry",
    function(x, type = 1L, sample_id = NULL) standardGeneric("annotGeometry")
)

#' @export
setGeneric(
    "annotGeometry<-",
    function(x, type = 1L, sample_id = NULL, translate = TRUE, ..., value)
        standardGeneric("annotGeometry<-")
)

#' @export
setGeneric(
    "annotGeometryNames",
    function(x) standardGeneric("annotGeometryNames")
)

#' @export
setGeneric(
    "annotGeometryNames<-",
    function(x, value) standardGeneric("annotGeometryNames<-")
)

#' @export
setGeneric("annotGeometries", function(x) standardGeneric("annotGeometries"))

#' @export
setGeneric(
    "annotGeometries<-",
    function(x, translate = TRUE, ..., value) {
        standardGeneric("annotGeometries<-")
    }
)

# Coercion--------
#' @export
setGeneric(
    "toSpatialFeatureExperiment",
    function(x, ...) standardGeneric("toSpatialFeatureExperiment")
)

# Spatial graphs--------
#' @export
setGeneric(
    "spatialGraphs",
    function(x, MARGIN = NULL, sample_id = "all", name = "all") {
        standardGeneric("spatialGraphs")
    }
)

#' @export
setGeneric(
    "spatialGraphs<-",
    function(x, MARGIN = NULL, sample_id = "all", name = "all", value) {
        standardGeneric("spatialGraphs<-")
    }
)

#' @export
setGeneric(
    "spatialGraph",
    function(x, type = 1L, MARGIN, sample_id = NULL) {
        standardGeneric("spatialGraph")
    }
)

#' @export
setGeneric(
    "spatialGraph<-",
    function(x, type = 1L, MARGIN, sample_id = NULL, value) {
        standardGeneric("spatialGraph<-")
    }
)

#' @export
setGeneric(
    "spatialGraphNames",
    function(x, MARGIN, sample_id = NULL) {
        standardGeneric("spatialGraphNames")
    }
)

#' @export
setGeneric(
    "spatialGraphNames<-",
    function(x, MARGIN, sample_id = NULL, value) {
        standardGeneric("spatialGraphNames<-")
    }
)

#' @export
setGeneric(
    "findSpatialNeighbors",
    function(x, ...) standardGeneric("findSpatialNeighbors")
)

# Have this generic and S4 dispatch to avoid conflict with sp::bbox
# Will change if this causes trouble
#' @export
setGeneric("bbox", function(sfe, sample_id = "all") standardGeneric("bbox"))

# Local results--------------
#' @export
setGeneric(
    "localResults",
    function(x, sample_id = "all", name = "all", features = NULL,
             colGeometryName = NULL, annotGeometryName = NULL,
             withDimnames = TRUE, swap_rownames = NULL, ...) {
        standardGeneric("localResults")
    }
)

#' @export
setGeneric(
    "localResults<-",
    function(x, sample_id = "all", name = "all", features = NULL,
             colGeometryName = NULL, annotGeometryName = NULL,
             withDimnames = TRUE, swap_rownames = NULL, ...,
             value) {
        standardGeneric("localResults<-")
    }
)

#' @export
setGeneric(
    "localResult",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, sample_id = NULL,
             withDimnames = TRUE, simplify = TRUE, swap_rownames = NULL) {
        standardGeneric("localResult")
    }
)

#' @export
setGeneric(
    "localResult<-",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, sample_id = NULL, swap_rownames = NULL,
             withDimnames = TRUE, value) {
        standardGeneric("localResult<-")
    }
)

#' @export
setGeneric("localResultNames", function(x) standardGeneric("localResultNames"))

#' @export
setGeneric(
    "localResultNames<-",
    function(x, value) standardGeneric("localResultNames<-")
)

#' @export
setGeneric(
    "localResultFeatures",
    function(x, type = 1L, colGeometryName = NULL, annotGeometryName = NULL,
             swap_rownames = NULL) {
        standardGeneric("localResultFeatures")
    }
)

#' @export
setGeneric(
    "localResultAttrs",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, swap_rownames = NULL) {
        standardGeneric("localResultAttrs")
    }
)

#' @export
setGeneric("unit", function(x) standardGeneric("unit"))

#' @export
setGeneric("transposeImg", function(x, ...) standardGeneric("transposeImg"))

if (!isGeneric("saveRDS")) {setGeneric("saveRDS", function (object, file="", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL) standardGeneric("saveRDS"))}
