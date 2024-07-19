#' Change sample IDs
#'
#' Change sample IDs in all fields of the SFE object where sample IDs are
#' present, not just the colData.
#'
#' @inheritParams sampleIDs
#' @param replacement A named character vector whose names are the existing
#'   sample IDs to be changed and whose values are the corresponding
#'   replacements.
#' @return An SFE object.
#' @concept Utilities
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
#' sampleIDs(sfe)
changeSampleIDs <- function(sfe, replacement) {
    for (i in seq_along(replacement)) {
        original <- names(replacement)[i]
        colData(sfe)$sample_id[colData(sfe)$sample_id == original] <-
            replacement[i]
        gs_names <- names(int_metadata(sfe)$spatialGraphs)
        names(int_metadata(sfe)$spatialGraphs)[gs_names == original] <-
            replacement[i]
        if (length(int_metadata(sfe)$annotGeometries)) {
            for (n in names(int_metadata(sfe)$annotGeometries)) {
                ag <- int_metadata(sfe)$annotGeometries[[n]]
                ind <- ag$sample_id == original
                int_metadata(sfe)$annotGeometries[[n]]$sample_id[ind] <- replacement[i]
            }
        }
        if (length(rowGeometries(sfe))) {
            nms <- rowGeometryNames(sfe)
            nms <- gsub(paste0(original, "$"), replacement[i], nms)
            rowGeometryNames(sfe) <- nms
            # Edge case: what if one sample_id includes another one?
            # e.g. sample01_x and x
        }
        if (nrow(imgData(sfe))) {
            imgData(sfe)$sample_id[imgData(sfe)$sample_id == original] <-
                replacement[i]
        }
        # Check spatial results
        # rowData
        rd_ind <- grepl(paste0("_", original), names(rowData(sfe)))
        if (any(rd_ind)) {
            names(rowData(sfe)) <- gsub(paste0("_", original), paste0("_", replacement[i]),
                                        names(rowData(sfe)))
        }
        # featureData
        if (!is.null(colFeatureData(sfe))) {
            nms <- names(colFeatureData(sfe))
            nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
            names(colFeatureData(sfe)) <- nms
        }
        # colGeometries
        for (n in colGeometryNames(sfe)) {
            if (!is.null(geometryFeatureData(sfe, n, 2L))) {
                nms <- names(geometryFeatureData(sfe, n, 2L))
                nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
                names(geometryFeatureData(sfe, n, 2L)) <- nms
            }
        }
        # annotGeometries
        for (n in annotGeometryNames(sfe)) {
            if (!is.null(geometryFeatureData(sfe, n, 3L))) {
                nms <- names(geometryFeatureData(sfe, n, 3L))
                nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
                names(geometryFeatureData(sfe, n, 3L)) <- nms
            }
        }
        # reducedDims
        for (n in reducedDimNames(sfe)) {
            if (!is.null(reducedDimFeatureData(sfe, n))) {
                nms <- names(reducedDimFeatureData(sfe, n))
                nms <- gsub(paste0("_", original), paste0("_", replacement[i]), nms)
                names(reducedDimFeatureData(sfe, n)) <- nms
            }
        }
    }
    sfe
}
# TODO: change sample ID based on geometry, which applies to all geometries.
# Need to note in documentation that featureData's results no longer apply.
