#' Get all unique sample IDs
#'
#' The title is self-explanatory.
#'
#' @param sfe A \code{SpatialFeatureExperiment} object.
#' @return A character vector of all unique entries of the \code{sample_id}
#' column in \code{colData(x)}.
#' @export
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sampleIDs(sfe)
sampleIDs <- function(sfe) unique(colData(sfe)$sample_id)

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
#' @export
#' @examples
#' library(SFEData)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
#' sampleIDs(sfe)
changeSampleIDs <- function(sfe, replacement) {
    for (i in seq_along(replacement)) {
        colData(sfe)$sample_id[colData(sfe)$sample_id == names(replacement)[i]] <-
            replacement[i]
        gs_names <- names(int_metadata(sfe)$spatialGraphs)
        names(int_metadata(sfe)$spatialGraphs)[gs_names == names(replacement)[i]] <-
            replacement[i]
        if (length(int_metadata(sfe)$annotGeometries)) {
            for (n in names(int_metadata(sfe)$annotGeometries)) {
                ag <- int_metadata(sfe)$annotGeometries[[n]]
                ind <- ag$sample_id == names(replacement)[i]
                int_metadata(sfe)$annotGeometries[[n]]$sample_id[ind] <- replacement[i]
            }
        }
    }
    sfe
}

.translate_value <- function(x, translate, value) {
    if (translate && !is.null(int_metadata(x)$orig_bbox)) {
        if (anyNA(value$sample_id) && nrow(value) == ncol(x))
            value$sample_id <- colData(x)$sample_id
        orig_bbox <- int_metadata(x)$orig_bbox
        # Don't translate if already translated
        curr_bbox <- st_bbox(value)
        samples <- unique(value$sample_id)
        if (length(samples) > 1L) {
            value$ID_ <- seq_len(nrow(value)) # Unlikely name
            df <- value[,c("ID_", "sample_id", "geometry")]
            df_split <- split(df, value$sample_id)
            df_split <- lapply(samples, function(s) {
                out <- df_split[[s]]
                og <- out$geometry - orig_bbox[c("xmin", "ymin"), s]
                bb <- st_bbox(og)
                if ((bb["xmin"] < 0 || bb["ymin"] < 0) || anyNA(bb)) return(out)
                else
                    out$geometry <- og
                out
            })
            df <- do.call(rbind, df_split)
            df <- df[match(value$ID_, df$ID_),]
            value$geometry <- df$geometry
            value$ID_ <- NULL
        } else {
            if (curr_bbox["xmin"] - orig_bbox["xmin", samples] < 0 ||
                curr_bbox["ymin"] - orig_bbox["ymin", samples] < 0) {
                return(value)
            }
            value$geometry <- value$geometry - orig_bbox[c("xmin", "ymin"), samples]
        }
    }
    value
}

#' Call metadata of an object
#'
#' @param object A \code{SpatialFeatureExperiment}, \code{SpatialExperiment}, or \code{SummarizedExperiment} object.
#' @return A \code{data.frame} object
#' @export
#' @importFrom SummarizedExperiment colData colData<- rowData
#' @examples
#' library(SFEData)
#' library(magrittr)
#' sfe <- McKellarMuscleData(dataset = "small")
#' sfe %>% callMeta() %>% str
callMeta <- function(object = NULL) {
    return(colData(object)@listData %>% as.data.frame.list())
}
