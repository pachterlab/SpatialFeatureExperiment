#' Save SpatialFeatureExperiment as RDS file
#'
#' Saving SFE objects as RDS files is complicated by the \code{SpatRaster} class
#' of the images. If present, the images need to be wrapped with the \code{\link{wrap}}
#' function in \code{terra} before serializing the SFE object. Otherwise the
#' images will be invalid pointers when the RDS is reloaded. If the image does
#' not fit in memory and its file source is unknown, then it will be written to
#' a temporary file, which is reloaded when the RDS file is loaded. When an SFE
#' object with images is read from an RDS file, the images will not be unwrapped
#' until necessary.
#'
#' @inheritParams base::saveRDS
#' @param object A \code{SpatialFeatureExperiment} object.
#' @return Invisibly \code{NULL}.
#' @importFrom terra wrap unwrap
#' @export
#' @examples
#' outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
#' samples <- file.path(outdir, paste0("sample0", 1:2))
#' sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
#' saveRDS(sfe, "foo.rds")
#'
setMethod("saveRDS", "SpatialFeatureExperiment",
          function(object, file = "", ascii = FALSE, version = NULL,
                   compress = TRUE, refhook = NULL) {
              if (is.null(imgData(object)))
                  base::saveRDS(object, file = file, ascii = ascii,
                                version = version, compress = compress,
                                refhook = refhook)
              else {
                  for (i in seq_len(nrow(imgData(object)))) {
                      img <- int_metadata(object)$imgData$data[[i]]
                      if (is(imgRaster(img), "SpatRaster"))
                          img@image <- wrap(img@image)
                      int_metadata(object)$imgData$data[[i]] <- img
                      # Then I need to add unwrap to functions that use the image
                  }
                  base::saveRDS(object, file = file, ascii = ascii,
                                version = version, compress = compress,
                                refhook = refhook)
              }
          })
