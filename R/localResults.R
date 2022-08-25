#' Organize results from local spatial statistics
#'
#' Local spatial statics like local Moran's I, local Geary's C, Getis-Ord Gi*,
#' and geographically weighted PCA loadings return values at each spatial
#' location. Just like dimension reductions, these results are clearly
#' associated with the braoder SFE object, so they should have a place within
#' the object. The \code{localResults} field in the SFE object stores these
#' results that has a value for each spatial location.
#'
#' @inheritParams dimGeometries
#' @name localResults
NULL

#' @rdname localResults
#' @export
setMethod("localResults", "SpatialFeatureExperiment",
          function(x, withDimnames = TRUE) {
            .get_intdimdata_all(x, MARGIN = 2L, withDimnames = withDimnames,
                                getfun = int_colData, key = "localResults")
          })

#' @rdname localResults
#' @export
setReplaceMethod("localResults", "SpatialFeatureExperiment",
                 function(x, withDimnames = TRUE, value) {
                   .set_intdimdata_all(x, MARGIN = 2L, withDimnames = withDimnames,
                                       translate = FALSE, sf = FALSE,
                                       getfun = int_colData,
                                       setfun = `int_colData<-`,
                                       key = "localResults",
                                       xdimfun = ncol,
                                       funstr = "localResults",
                                       xdimstr = "ncol", value)
                 })

#' @rdname localResults
#' @export
setMethod("localResultNames", "SpatialFeatureExperiment",
          function(x) {
            .get_internal_names(x,
                                getfun=int_colData,
                                key="localResults")
          })

#' @rdname localResults
#' @export
setReplaceMethod("localResultNames",
                 c("SpatialFeatureExperiment", "character"),
                 function(x, value) {
                   .set_internal_names(x, value,
                                       getfun=int_colData,
                                       setfun=`int_colData<-`,
                                       key="localResults")
                 })

#' @rdname localResults
#' @export
setMethod("localResult", c("SpatialFeatureExperiment", "missing"),
          function(x, type, sample_id = NULL, withDimnames = TRUE) {
            .get_internal_missing(x,
                                  basefun=localResult,
                                  namefun=localResultNames,
                                  funstr="localResult",
                                  withDimnames=withDimnames,
                                  sample_id = sample_id)
            })

#' @rdname localResults
#' @export
setMethod("localResult", c("SpatialFeatureExperiment", "numeric"),
          function(x, type, sample_id = NULL, withDimnames = TRUE) {
            .get_internal_id(x, type, MARGIN = 2L, sample_id = sample_id,
                             withDimnames = withDimnames,
                             .get_internal_fun = .get_internal_integer,
                             getfun = int_colData,
                             key = "localResults", funstr = "localResult",
                             substr = "type")
          })

#' @rdname localResults
#' @export
setMethod("localResult", c("SpatialFeatureExperiment", "character"),
          function(x, type, sample_id = NULL, withDimnames = TRUE) {
            .get_internal_id(x, type, MARGIN = 2L, sample_id = sample_id,
                             withDimnames = withDimnames,
                             .get_internal_fun = .get_internal_character,
                             getfun = int_colData,
                             key = "localResults", funstr = "localResult",
                             substr = "type")
          })

#' @rdname localResults
#' @export
setReplaceMethod("localResult", c("SpatialFeatureExperiment", "missing"),
                 function(x, type, sample_id = NULL, withDimnames=TRUE, value) {
                   .set_internal_missing(x, value,
                                         withDimnames=withDimnames,
                                         sample_id = sample_id,
                                         basefun=`localResult<-`,
                                         namefun=localResultNames)
                 })

#' @rdname localResults
#' @export
setReplaceMethod("localResult", c("SpatialFeatureExperiment", "numeric"),
                 function(x, type, sample_id = NULL, withDimnames=TRUE, value) {
                   .set_internal_id(x, type, MARGIN = 2L, sample_id = sample_id,
                                    withDimnames = withDimnames,
                                    translate = FALSE, sf = FALSE,
                                    .set_internal_fun = .set_internal_numeric,
                                    getfun = int_colData,
                                    setfun = `int_colData<-`,
                                    key = "localResults",
                                    xdimfun = ncol,
                                    funstr = "localResult", xdimstr = "ncol",
                                    substr = "type", value)
                 })

#' @rdname localResults
#' @export
setReplaceMethod("localResult", c("SpatialFeatureExperiment", "character"),
                 function(x, type, sample_id = NULL, withDimnames=TRUE, value) {
                   .set_internal_id(x, type, MARGIN = 2L, sample_id = sample_id,
                                    withDimnames = withDimnames,
                                    translate = FALSE, sf = FALSE,
                                    .set_internal_fun = .set_internal_character,
                                    getfun = int_colData,
                                    setfun = `int_colData<-`,
                                    key = "localResults",
                                    xdimfun = ncol,
                                    funstr = "localResult", xdimstr = "ncol",
                                    substr = "type", value)
                 })
