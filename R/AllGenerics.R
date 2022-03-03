# Geometry getters and setters should be kind of similar to implementation of
# reducedDims in SCE

# Getters and setters for dimGeometry--------
# For colGeometry and rowGeometry
#' @export
setGeneric("dimGeometry", function(x, type, MARGIN, withDimnames=TRUE) standardGeneric("dimGeometry"))

#' @export
setGeneric("dimGeometry<-", function(x, type, MARGIN, withDimnames=TRUE,..., value) standardGeneric("dimGeometry<-"))

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
setGeneric("annotGeometry", function(x, type) standardGeneric("annotGeometry"))

#' @export
setGeneric("annotGeometry<-", function(x, type, ..., value) standardGeneric("annotGeometry<-"))

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
setGeneric("spatialGraphs", function(x, MARGIN) standardGeneric("spatialGraphs"))

#' @export
setGeneric("spatialGraphs<-", function(x, MARGIN, value) standardGeneric("spatialGraphs<-"))

#' @export
setGeneric("spatialGraph", function(x, type, MARGIN) standardGeneric("spatialGraph"))

#' @export
setGeneric("spatialGraph<-", function(x, type, MARGIN, value) standardGeneric("spatialGraph<-"))

# Operations that only apply to SFE due to geometries------
#' @export
setGeneric("crop", function(x, y) standardGeneric("crop"))
