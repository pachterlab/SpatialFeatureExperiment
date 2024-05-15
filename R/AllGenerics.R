# Geometry getters and setters should be kind of similar to implementation of
# reducedDims in SCE

# Getters and setters for dimGeometry--------
# For colGeometry and rowGeometry
#' @export
setGeneric(
    "dimGeometry",
    function(x, type = 1L, MARGIN, sample_id = 1L,
             withDimnames = TRUE) {
        standardGeneric("dimGeometry")
    }
)

#' @export
setGeneric("dimGeometry<-", function(x, type = 1L, MARGIN, sample_id = 1L,
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
    function(x, type = 1L, sample_id = 1L) standardGeneric("annotGeometry")
)

#' @export
setGeneric(
    "annotGeometry<-",
    function(x, type = 1L, sample_id = 1L, translate = TRUE, ..., value)
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
    function(x, type = 1L, MARGIN, sample_id = 1L) {
        standardGeneric("spatialGraph")
    }
)

#' @export
setGeneric(
    "spatialGraph<-",
    function(x, type = 1L, MARGIN, sample_id = 1L, value) {
        standardGeneric("spatialGraph<-")
    }
)

#' @export
setGeneric(
    "spatialGraphNames",
    function(x, MARGIN, sample_id = 1L) {
        standardGeneric("spatialGraphNames")
    }
)

#' @export
setGeneric(
    "spatialGraphNames<-",
    function(x, MARGIN, sample_id = 1L, value) {
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
setGeneric("bbox", function(sfe, sample_id = "all", ...) standardGeneric("bbox"))

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
             annotGeometryName = NULL, sample_id = 1L,
             withDimnames = TRUE, simplify = TRUE, swap_rownames = NULL) {
        standardGeneric("localResult")
    }
)

#' @export
setGeneric(
    "localResult<-",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, sample_id = 1L, swap_rownames = NULL,
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

setGeneric(".transpose_img", function(x, bbox_all, ...) standardGeneric(".transpose_img"))

setGeneric(".mirror_img", function(x, direction, bbox_all, ...) standardGeneric(".mirror_img"))

setGeneric(".rotate_img", function(x, degrees, bbox_all, ...) standardGeneric(".rotate_img"))

#' @export
setGeneric("scaleImg", function(x, factor, ...) standardGeneric("scaleImg"))

#' @export
setGeneric("translateImg", function(x, v, ...) standardGeneric("translateImg"))

#' @export
setGeneric("affineImg", function(x, M, v, ...) standardGeneric("affineImg"))

#' @export
setGeneric("cropImg", function(x, bbox, ...) standardGeneric("cropImg"))

if (!isGeneric("saveRDS")) {setGeneric("saveRDS", function (object, file="", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL) standardGeneric("saveRDS"))}

#' @export
setGeneric("toExtImage", function(x, ...) standardGeneric("toExtImage"))

#' @export
setGeneric("toSpatRasterImage", function(x, ...) standardGeneric("toSpatRasterImage"))

#' @export
setGeneric("isFull", function(x) standardGeneric("isFull"))

#' @export
setGeneric("origin", function(x) standardGeneric("origin"))

setGeneric("origin<-", function(x, value) standardGeneric("origin<-"))

#' @export
setGeneric("transformation", function(x) standardGeneric("transformation"))

setGeneric("transformation<-", function(x, value) standardGeneric("transformation<-"))

if (!isGeneric("unwrap")) {setGeneric("unwrap", function(x, ...) standardGeneric("unwrap"))}

#' @export
setGeneric("Img<-", function(x, sample_id = 1L, image_id, scale_fct = 1, value) standardGeneric("Img<-"))
