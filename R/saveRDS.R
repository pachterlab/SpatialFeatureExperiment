#' Save SpatialFeatureExperiment as RDS file
#'
#' Saving SFE objects as RDS files is complicated by the \code{SpatRaster} class
#' of the images. If present, the images need to be wrapped with the \code{\link[terra]{wrap}}
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
#' @importFrom BiocGenerics saveRDS
#' @export
#' @concept Utilities
#' @examples
#' outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
#' samples <- file.path(outdir, paste0("sample0", 1:2))
#' sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
#' saveRDS(sfe, "foo.rds")
#' # Clean up
#' file.remove("foo.rds")
setMethod("saveRDS", "SpatialFeatureExperiment",
          function(object, file = "", ascii = FALSE, version = NULL,
                   compress = TRUE, refhook = NULL) {
              if (containsOutOfMemoryData(object)) {
                  warning("This object contains out of memory data that may
                          break when the RDS file is moved.")
              }
              if (!nrow(imgData(object)))
                  base::saveRDS(object, file = file, ascii = ascii,
                                version = version, compress = compress,
                                refhook = refhook)
              else {
                  for (i in seq_len(nrow(imgData(object)))) {
                      img <- int_metadata(object)$imgData$data[[i]]
                      if (inherits(img, "SpatRasterImage"))
                          img <- new("PackedRasterImage", wrap(img))
                      int_metadata(object)$imgData$data[[i]] <- img
                  }
                  base::saveRDS(object, file = file, ascii = ascii,
                                version = version, compress = compress,
                                refhook = refhook)
              }
          })
# From terra
setMethod("readRDS", signature(file="character"),
          function (file = "", refhook = NULL) {
              x <- base::readRDS(file=file, refhook=refhook)
              unwrap(x)
          }
)

setMethod("unwrap", "SpatialFeatureExperiment",
          function(x) {
              for (i in seq_len(nrow(imgData(x)))) {
                  img <- int_metadata(x)$imgData$data[[i]]
                  if (inherits(img, "PackedRasterImage"))
                      img <- SpatRasterImage(unwrap(img))
                  else if (inherits(img, "SpatRasterImage")) {
                      old_slot <- tryCatch(img@image, error = function(e) NULL)
                      if (!is.null(old_slot)) {
                          if (inherits(old_slot, "SpatRaster")) img <- old_slot
                          if (inherits(old_slot, "PackedSpatRaster")) img <- unwrap(old_slot)
                          img <- SpatRasterImage(img)
                      }
                  }
                  int_metadata(x)$imgData$data[[i]] <- img
              }
              x
          })

#' Whether an SFE object contains out of memory data
#' 
#' Out of memory data, such as \code{DelayedArray}, some \code{SpatRasterImage}
#' objects, and \code{BioFormatsImage}, will break if saved as RDS. This method
#' of \code{\link[BiocGenerics]{containsOutOfMemoryData}} checks if an SFE 
#' object has out of memory data, specifically the images. Having out of memory
#' data will result into an error when \code{saveRDS} is called; we recommend
#' using the \code{alabaster.sfe} package instead.
#' 
#' @param object An SFE object
#' @return TRUE or FALSE
#' @importFrom BiocGenerics containsOutOfMemoryData
#' @export
#' @concept Utilities
#' @examples
#' outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
#' samples <- file.path(outdir, paste0("sample0", 1:2))
#' sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
#' containsOutOfMemoryData(sfe)
setMethod("containsOutOfMemoryData", "SpatialFeatureExperiment",
          function(object) {
              imgs <- imgData(object)
              imgData(object) <- NULL
              object <- as(object, "SingleCellExperiment")
              out_se <- containsOutOfMemoryData(object)
              is_img_om <- vapply(imgs$data, function(x) {
                  !is.na(imgSource(x))
              }, FUN.VALUE = logical(1))
              any(c(out_se, is_img_om))
          })

