# Geometry getters and setters should be kind of similar to implementation of
# reducedDims in SCE

# Getters and setters for dimGeometry--------
# For colGeometry and rowGeometry
#' @export
setGeneric("dimGeometry", function(x, type, MARGIN, sample_id, withDimnames=TRUE) standardGeneric("dimGeometry"))

#' @export
setGeneric("dimGeometry<-", function(x, type, MARGIN, sample_id, withDimnames=TRUE,..., value) standardGeneric("dimGeometry<-"))

#' @export
setGeneric("dimGeometryNames", function(x, MARGIN) standardGeneric("dimGeometryNames"))

#' @export
setGeneric("dimGeometryNames<-", function(x, MARGIN, value) standardGeneric("dimGeometryNames<-"))

#' @export
setGeneric("dimGeometries", function(x, MARGIN, withDimnames=TRUE) standardGeneric("dimGeometries"))

#' @export
setGeneric("dimGeometries<-", function(x, MARGIN, withDimnames=TRUE,..., value) standardGeneric("dimGeometries<-"))

# Getters and setters for annotGeometry---------
#' @export
setGeneric("annotGeometry", function(x, type, sample_id) standardGeneric("annotGeometry"))

#' @export
setGeneric("annotGeometry<-", function(x, type, sample_id, ..., value) standardGeneric("annotGeometry<-"))

#' @export
setGeneric("annotGeometryNames", function(x) standardGeneric("annotGeometryNames"))

#' @export
setGeneric("annotGeometryNames<-", function(x, value) standardGeneric("annotGeometryNames<-"))

#' @export
setGeneric("annotGeometries", function(x) standardGeneric("annotGeometries"))

#' @export
setGeneric("annotGeometries<-", function(x, ..., value) standardGeneric("annotGeometries<-"))

# Coercion--------
#' @export
setGeneric("toSpatialFeatureExperiment", function(x, ...) standardGeneric("toSpatialFeatureExperiment"))

# Spatial graphs--------
#' @export
setGeneric("spatialGraphs", function(x, MARGIN, sample_id) standardGeneric("spatialGraphs"))

#' @export
setGeneric("spatialGraphs<-", function(x, MARGIN, sample_id, value) standardGeneric("spatialGraphs<-"))

#' @export
setGeneric("colGraphs", function(x, sample_id) standardGeneric("colGraphs"))

#' @export
setGeneric("colGraphs<-", function(x, sample_id, value) standardGeneric("colGraphs<-"))

#' @export
setGeneric("rowGraphs", function(x, sample_id) standardGeneric("rowGraphs"))

#' @export
setGeneric("rowGraphs<-", function(x, sample_id, value) standardGeneric("rowGraphs<-"))

#' @export
setGeneric("annotGraphs", function(x, sample_id) standardGeneric("annotGraphs"))

#' @export
setGeneric("annotGraphs<-", function(x, sample_id, value) standardGeneric("annotGraphs<-"))

#' @export
setGeneric("spatialGraph", function(x, type, MARGIN, sample_id) standardGeneric("spatialGraph"))

#' @export
setGeneric("spatialGraph<-", function(x, type, MARGIN, sample_id, value) standardGeneric("spatialGraph<-"))

#' @export
setGeneric("spatialGraphNames", function(x, MARGIN, sample_id) standardGeneric("spatialGraphNames"))

#' @export
setGeneric("spatialGraphNames<-", function(x, MARGIN, sample_id, value) standardGeneric("spatialGraphNames<-"))

# Operations that only apply to SFE due to geometries------
#' @export
setGeneric("crop", function(x, y) standardGeneric("crop"))
