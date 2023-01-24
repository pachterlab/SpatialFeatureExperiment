#' Read 10X Visium data as SpatialFeatureExperiment
#'
#' Read Space Ranger output as a SpatialFeatureExperiment object, where spots
#' are represented with polygons in the colGeometry called "spotPoly". Other
#' geometries can be added later after the dataset is read. If \code{data =
#' "filtered"}, then spatial neighborhood graphs of the spots are also computed
#' and stored in the colGraph called "visium" in all samples for downstream
#' spatial analyses.
#'
#' @inheritParams SpatialExperiment::read10xVisium
#' @inheritParams findVisiumGraph
#' @inheritParams SpatialFeatureExperiment
#' @param type Either "HDF5", and the matrix will be represented as
#'   \code{TENxMatrix}, or "sparse", and the matrix will be read as a
#'   \code{dgCMatrix}.
#' @param dirs Directory for each sample that contains the \code{spatial} and
#'   \code{raw/filtered_featues_bc_matrix} directories. By default, the
#'   \code{outs} directory under the directory specified in the \code{samples}
#'   argument, as in Space Ranger output. Change the \code{dirs} argument if you
#'   have moved or renamed the output directory.
#' @note The \code{as(<dgTMatrix>, "dgCMatrix") is deprecated} warning comes
#'   from the \code{DropletUtils} package which is used by
#'   \code{SpatialExperiment} to read 10X outputs. This will be fixed when
#'   \code{SpatialExperiment} switches to TENxIO.
#' @importFrom SpatialExperiment read10xVisium
#' @importFrom rjson fromJSON
#' @importFrom SummarizedExperiment rowData<-
#' @return A SpatialFeatureExperiment object
#' @export
#' @examples
#' dir <- system.file("extdata", package = "SpatialFeatureExperiment")
#'
#' sample_ids <- c("sample01", "sample02")
#' samples <- file.path(dir, sample_ids, "outs")
#'
#' list.files(samples[1])
#' list.files(file.path(samples[1], "spatial"))
#' (sfe <- read10xVisiumSFE(samples, sample_ids,
#'     type = "sparse", data = "filtered",
#'     load = FALSE
#' ))
read10xVisiumSFE <- function(samples = "",
                             dirs = file.path(samples, "outs"),
                             sample_id = paste0(
                                 "sample",
                                 sprintf(
                                     "%02d",
                                     seq_along(samples)
                                 )
                             ),
                             type = c("HDF5", "sparse"),
                             data = c("filtered", "raw"),
                             images = "lowres",
                             load = TRUE, style = "W", zero.policy = NULL,
                             BPPARAM = SerialParam()) {
    # Read one sample at a time, in order to get spot diameter one sample at a time
    sfes <- lapply(seq_along(samples), function(i) {
        suppressWarnings(o <- read10xVisium(dirs[i], sample_id[i], type, data, images, load))
        # Transpose coordinates
        spatialCoords(o) <- spatialCoords(o)[, c(2, 1)]
        scalefactors <- fromJSON(file = file.path(
            dirs[i], "spatial",
            "scalefactors_json.json"
        ))
        o <- .spe_to_sfe(o,
            colGeometries = NULL, rowGeometries = NULL,
            annotGeometries = NULL, spatialCoordsNames = NULL,
            annotGeometryType = NULL, spatialGraphs = NULL,
            spotDiameter = scalefactors[["spot_diameter_fullres"]],
            unit = "full_res_image_pixels", BPPARAM = BPPARAM
        )
        # Add spatial enrichment if present
        fn <- file.path(dirs[i], "spatial", "spatial_enrichment.csv")
        if (file.exists(fn)) {
            enrichment <- read.csv(fn)
            row_inds <- match(rownames(o), enrichment$Feature.ID)
            # Let me not worry about different samples having different genes for now
            if (i == 1L) {
                rowData(o) <- cbind(rowData(o), 
                                    Feature.Type = enrichment[row_inds, "Feature.Type"])
            }
            cols_use <- c("I", "P.value", "Adjusted.p.value", 
                          "Feature.Counts.in.Spots.Under.Tissue", 
                          "Median.Normalized.Average.Counts", 
                          "Barcodes.Detected.per.Feature")
            enrichment2 <- enrichment[row_inds, cols_use]
            names(enrichment2) <- paste(names(enrichment2), sample_id[i], 
                                        sep = "_")
            rowData(o) <- cbind(rowData(o), enrichment2)
        }
        # Add barcode fluorescence intensity if present
        fn2 <- file.path(dirs[i], "spatial", "barcode_fluorescence_intensity.csv")
        if (file.exists(fn)) {
            fluo <- read.csv(fn2)
            row_inds <- match(colnames(o), fluo$barcode)
            fluo$barcode <- NULL
            fluo$in_tissue <- NULL
            colData(o) <- cbind(colData(o), fluo[row_inds,])
        }
        o
    })
    out <- do.call(cbind, sfes)
    if (data == "filtered") {
        colGraphs(out, sample_id = "all", name = "visium") <-
            findVisiumGraph(out,
                sample_id = "all", style = style,
                zero.policy = zero.policy
            )
    }
    out
}
