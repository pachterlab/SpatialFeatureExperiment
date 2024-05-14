#' Get and set results from local spatial statistics
#'
#' Local spatial statics like local Moran's I, local Geary's C, Getis-Ord Gi*,
#' and geographically weighted summary statistics return values at each spatial
#' location. Just like dimension reductions, these results are clearly
#' associated with the broader SFE object, so they should have a place within
#' the object. However, a separate field is needed because these analyses are
#' conceptually distinct from dimension reduction. Also, each feature (e.g.
#' gene) can have its own results with values at each location. The
#' \code{localResults} field in the SFE object stores these results that has a
#' value for each spatial location.
#'
#' @inheritParams dimGeometries
#' @inheritParams spatialGraphs
#' @param ... Ignored
#' @param type Name or index of the spatial method used, such as "localmoran".
#' @param name Name of the spatial method used, such as "localmoran".
#' @param value Values to set, should be either a matrix or a data frame.
#' @param feature Feature whose local results to get or set, for
#'   \code{localResult} getter and setter for one feature at a time.
#' @param features Features whose local results to get or set, for
#'   \code{localResults} getter and setter for multiple features at a time.
#' @param colGeometryName Which \code{colGeometry} to get or set local results.
#' @param annotGeometryName Which \code{annotGeometry} to get or set local
#'   results.
#' @param simplify Basically whether to return the content of the list rather
#'   than a list when the list only has one element, such as results for one
#'   type and one feature.
#' @param swap_rownames Name of a column in \code{rowData} to identify features
#'   instead of the row names of the SFE object. For example, if the row names
#'   of the SFE object are Ensembl IDs and gene symbols are in the "symbol"
#'   column in \code{rowData}, then putting "symbol" for this argument will use
#'   the gene symbols to identify which gene's local results to get or set.
#' @aliases localResults localResults<- localResult localResult<-
#'   localResultNames localResultNames<- localResultFeatures localResultAttrs
#' @return \code{localResults} returns a named list each element of which is a
#'   set of local results of interest. \code{localResult} returns a matrix or a
#'   data frame, whichever the original is when it's set.
#'   \code{localResultNames} returns a character vector. Setters return an SFE
#'   object with the desired field set. For genes and \code{colData} columns,
#'   the local results are stored in the \code{localResults} field in
#'   \code{int_colData}, whereas for \code{colGeometries} and
#'   \code{annotGeometries}, the local results are stored as columns in the same
#'   \code{sf} data frames. \code{localResultFeatures} returns a character
#'   vector of names of features for which local results are available.
#'   \code{localResultAttrs} returns a character vector of the column names of
#'   the local results of one type for one feature. It returns \code{NULL} if
#'   the results are a vector.
#' @docType methods
#' @name localResults
#' @concept Getters and setters
#' @examples
#' # Toy example
#' sfe <- readRDS(system.file("extdata/sfe_toy.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' # localResults functions are written for organizing results from local
#' # spatial statistics (see the Voyager package). But for the examples here,
#' # random toy matrices are used. The real results are often matrices, with a
#' # matrix for each feature.
#' library(S4Vectors)
#' set.seed(29)
#' toy_res1 <- matrix(rnorm(10),
#'     nrow = 5, ncol = 2,
#'     dimnames = list(colnames(sfe), c("meow", "purr"))
#' )
#' toy_res1b <- matrix(rgamma(10, shape = 2),
#'     nrow = 5, ncol = 2,
#'     dimnames = list(colnames(sfe), c("meow", "purr"))
#' )
#' toy_df1 <- DataFrame(gene1 = I(toy_res1), gene2 = I(toy_res1b))
#'
#' toy_res2 <- matrix(rpois(10, lambda = 2),
#'     nrow = 5, ncol = 2,
#'     dimnames = list(colnames(sfe), c("sassy", "tortitude"))
#' )
#' toy_df2 <- DataFrame(gene1 = I(toy_res2))
#' # Set all local results
#' localResults(sfe) <- list(localmoran = toy_df1, Gistar = toy_df2)
#' # Get all local results
#' lrs <- localResults(sfe)
#'
#' # Set results of the same type for multiple genes
#' localResults(sfe, name = "localmoran") <- toy_df1
#' # Can also use a list
#' localResults(sfe, name = "localmoran") <- as.list(toy_df1)
#' # Get results of the same type for multiple genes
#' lrs <- localResults(sfe, name = "localmoran", features = c("gene1", "gene2"))
#'
#' # Set results for one type and one gene
#' localResult(sfe, "localmoran", feature = "gene1") <- toy_res1
#' # Get results for one type and one gene
#' lr <- localResult(sfe, "localmoran", feature = "gene1")
#'
#' # Set results for a feature in colGeometries
#' cg_toy <- readRDS(system.file("extdata/cg_toy.rds",
#'     package = "SpatialFeatureExperiment"
#' ))
#' colGeometry(sfe, "cg") <- cg_toy
#' localResult(sfe, "localmoran",
#'     feature = "gene1",
#'     colGeometryName = "cg"
#' ) <- toy_res1
#' # Get results for a feature in colGeometries
#' lr <- localResult(sfe, "localmoran", "gene1", colGeometryName = "cg")
NULL

# Get all results for all features and all types
#' @rdname localResults
#' @export
setMethod("localResults", "SpatialFeatureExperiment",
          function(x, sample_id = "all", name = "all", features = NULL,
                   colGeometryName = NULL,
                   annotGeometryName = NULL, withDimnames = TRUE,
                   swap_rownames = NULL, ...) {
              if (isTRUE(all.equal(sample_id, "all")) && name == "all") {
                  .get_intdimdata_all(x,
                                      MARGIN = 2L, withDimnames = withDimnames,
                                      getfun = int_colData, key = "localResults"
                  )
              } else {
                  localResult(x, name,
                              feature = features, sample_id = sample_id,
                              withDimnames = withDimnames,
                              colGeometryName = colGeometryName,
                              annotGeometryName = annotGeometryName,
                              simplify = FALSE, swap_rownames = swap_rownames
                  )
              }
          })

#' @rdname localResults
#' @export
setReplaceMethod("localResults", "SpatialFeatureExperiment",
                 function(x, sample_id = "all", name = "all", features = NULL,
                          colGeometryName = NULL, annotGeometryName = NULL,
                          withDimnames = TRUE, swap_rownames = NULL,..., value) {
                     if (isTRUE(all.equal(sample_id, "all")) && name == "all") {
                         .set_intdimdata_all(x,
                                             MARGIN = 2L,
                                             withDimnames = withDimnames,
                                             translate = FALSE, sf = FALSE,
                                             getfun = int_colData,
                                             setfun = `int_colData<-`,
                                             key = "localResults",
                                             xdimfun = ncol,
                                             funstr = "localResults",
                                             xdimstr = "ncol", value
                         )
                     } else {
                         `localResult<-`(x, name,
                                         feature = features,
                                         sample_id = sample_id,
                                         withDimnames = withDimnames,
                                         colGeometryName = colGeometryName,
                                         annotGeometryName = annotGeometryName,
                                         value = value
                         )
                     }
                 })

# Which other scenario to get or set multiple results?
# When computing or plotting results of the same metric for multiple genes

#' @rdname localResults
#' @export
setMethod(
    "localResultNames", "SpatialFeatureExperiment",
    function(x) {
        .get_internal_names(x,
            getfun = int_colData,
            key = "localResults"
        )
    }
)

#' @rdname localResults
#' @export
setReplaceMethod(
    "localResultNames",
    c("SpatialFeatureExperiment", "character"),
    function(x, value) {
        .set_internal_names(x, value,
            getfun = int_colData,
            setfun = `int_colData<-`,
            key = "localResults"
        )
    }
)

#' @rdname localResults
#' @export
setMethod(
    "localResultFeatures", "SpatialFeatureExperiment",
    function(x, type = 1L, colGeometryName = NULL,
             annotGeometryName = NULL, swap_rownames = NULL) {
        tryCatch(names(localResults(x, "all", type,
            colGeometryName = colGeometryName,
            annotGeometryName = annotGeometryName, withDimnames = FALSE,
            swap_rownames = swap_rownames
        )), error = function(e) NULL)
    }
)

#' @rdname localResults
#' @export
setMethod(
    "localResultAttrs", "SpatialFeatureExperiment",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, swap_rownames = NULL) {
        tryCatch(colnames(localResult(x, type, feature, colGeometryName,
            annotGeometryName, sample_id = "all", withDimnames = FALSE,
            swap_rownames = swap_rownames
        )), error = function(e) NULL)
    }
)

# Here "feature" can be a character vector for multiple features, but it will
# not be documented. Use localResults for that.
#' @rdname localResults
#' @export
setMethod(
    "localResult", "SpatialFeatureExperiment",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, sample_id = 1L,
             withDimnames = TRUE, simplify = TRUE, swap_rownames = NULL) {
        .get_internal_feature(x,
            type = type, feature = feature,
            MARGIN = 2L,
            colGeometryName = colGeometryName,
            annotGeometryName = annotGeometryName,
            sample_id = sample_id,
            withDimnames = withDimnames,
            .get_internal_fun = .get_internal,
            getfun = int_colData,
            key = "localResults", funstr = "localResult",
            substr = "type", namestr = "localResult", simplify = simplify,
            swap_rownames = swap_rownames
        )
    }
)

#' @rdname localResults
#' @export
setReplaceMethod(
    "localResult", "SpatialFeatureExperiment",
    function(x, type = 1L, feature, colGeometryName = NULL,
             annotGeometryName = NULL, sample_id = 1L,
             withDimnames = TRUE, value) {
        .set_internal_feature(x, type,
            feature = feature,
            colGeometryName = colGeometryName,
            annotGeometryName = annotGeometryName,
            MARGIN = 2L, sample_id = sample_id,
            withDimnames = withDimnames,
            translate = FALSE, sf = FALSE,
            .get_all_fun = localResults,
            .set_all_fun = `localResults<-`,
            .set_internal_fun = .set_internal,
            getfun = int_colData,
            setfun = `int_colData<-`,
            key = "localResults",
            xdimfun = ncol,
            funstr = "localResult",
            xdimstr = "ncol",
            substr = "type", value = value
        )
    }
)
