.concat_annotgeometries <- function(x, y) {
    # Check the names of the geometries in each object
    # For geometries with the same name, check column names for rbind
    # If the columns are different, then append sample_id to the geometry names
    ag_names_x <- annotGeometryNames(x)
    ag_names_y <- annotGeometryNames(y)
    if (is.null(x) && is.null(y)) {
        return(NULL)
    }
    names_intersect <- intersect(ag_names_x, ag_names_y)
    if (length(names_intersect)) {
        # A bit of repetition, just to get the names right
        which_rbind <- vapply(names_intersect, function(n) {
            ag_x <- annotGeometry(x, n)
            ag_y <- annotGeometry(y, n)
            setequal(names(ag_x), names(ag_y))
        }, FUN.VALUE = logical(1))

        out_rbind <- lapply(names_intersect[which_rbind], function(n) {
            ag_x <- annotGeometry(x, n)
            ag_y <- annotGeometry(y, n)
            ag_y <- ag_y[, names(ag_x)]
            rbind(ag_x, ag_y)
        })
        out_append <- lapply(names_intersect[!which_rbind], function(n) {
            ag_x <- annotGeometry(x, n)
            ag_y <- annotGeometry(y, n)
            setNames(list(ag_x, ag_y), paste(n,
                c(ag_x$sample_id[1], ag_y$sample_id[1]),
                sep = "_"
            ))
        })
        out_append <- unlist(out_append, recursive = FALSE)
        names(out_rbind) <- names_intersect[which_rbind]
        out_x <- annotGeometries(x)[setdiff(ag_names_x, ag_names_y)]
        out_y <- annotGeometries(y)[setdiff(ag_names_y, ag_names_x)]
    } else {
        out_rbind <- out_append <- NULL
        out_x <- annotGeometries(x)
        out_y <- annotGeometries(y)
    }
    out <- c(out_rbind, out_append, out_x, out_y)
    return(out)
}

#' Concatenate SpatialFeatureExperiment objects
#'
#' On top of the \code{cbind} method of \code{SpatialExperiment}, this method is
#' needed to properly merge the \code{spatialGraphs} field in the different SFE
#' objects.
#'
#' @param ... SFE objects to cbind.
#' @param deparse.level See \code{?\link[base]{rbind}}.
#' @return A combined SFE object.
#' @importFrom BiocGenerics cbind
#' @export
#' @examples
#' library(SFEData)
#' sfe_small <- McKellarMuscleData(dataset = "small")
#' sfe_small2 <- McKellarMuscleData(dataset = "small2")
#' sfe2 <- cbind(sfe_small, sfe_small2)
setMethod(
    "cbind", "SpatialFeatureExperiment",
    function(..., deparse.level = 1) {
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }
        args <- list(...)
        out <- do.call(
            callNextMethod,
            c(args, list(deparse.level = deparse.level))
        )
        colnames(out) <- make.unique(colnames(out), sep = "-")
        # Combine the annotGeometries
        has_ag <- vapply(args, function(a) length(annotGeometries(a)) > 0L,
            FUN.VALUE = logical(1)
        )
        if (any(has_ag)) {
            new_ag <- Reduce(.concat_annotgeometries, args)
            int_metadata(out)[["annotGeometries"]] <- NULL
            int_metadata(out)[["annotGeometries"]] <- new_ag
        }
        # Combine the spatialGraphs
        has_sg <- vapply(args, function(a) !is.null(spatialGraphs(a)),
            FUN.VALUE = logical(1)
        )
        if (any(has_sg)) {
            sgs_use <- lapply(args, function(a) {
                if (is.null(spatialGraphs(a))) {
                    .initialize_spatialGraphs(a)
                } else {
                    int_metadata(a)[["spatialGraphs"]]
                }
            })
            new_sgs <- do.call(cbind, sgs_use)
            int_metadata(out)[["spatialGraphs"]] <- NULL
            int_metadata(out)$spatialGraphs <- new_sgs
        }
        # Clean up duplicate version field in int_metadata
        # Bug in SCE
        im_names <- unique(names(int_metadata(out)))
        names(int_metadata(out)) <- make.unique(names(int_metadata(out)))
        int_metadata(out) <- int_metadata(out)[im_names]
        out
    }
)
