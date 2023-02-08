library(SpatialExperiment)
library(SingleCellExperiment)
library(DropletUtils)
library(sf)

dir <- system.file("extdata/sample01", package = "SpatialFeatureExperiment")

sce <- read10xCounts(file.path(dir, "outs", "filtered_feature_bc_matrix"))

# Read in spatial coordinates
coords <- read.csv(file.path(dir, "outs", "spatial", "tissue_positions.csv"))
int_colData(sce)$spatialCoords <- as.matrix(coords[,c("pxl_row_in_fullres", "pxl_col_in_fullres")])
spe <- as(sce, "SpatialExperiment")

test_that("Convert SPE and SCE to SFE", {
    sfe1 <- as(spe, "SpatialFeatureExperiment")
    sfe2 <- as(sce, "SpatialFeatureExperiment")
    sfe3 <- toSpatialFeatureExperiment(spe)
    colData(sce) <- cbind(colData(sce), coords[,c("pxl_row_in_fullres", "pxl_col_in_fullres")])
    sfe4 <- toSpatialFeatureExperiment(sce, spatialCoordsNames = c("pxl_row_in_fullres",
                                                                   "pxl_col_in_fullres"))
    sfe <- read10xVisiumSFE(dir, type = "sparse", data = "filtered")
    centroids_check <- st_centroid(spotPoly(sfe))

    expect_s4_class(sfe1, "SpatialFeatureExperiment")
    expect_equal(sfe2, sfe1)
    expect_equal(sfe3, sfe1)
    expect_equal(sfe4, sfe1)
    expect_equal(centroids(sfe1), centroids_check, ignore_attr = "row.names")
    expect_equal(spatialCoords(sfe1), spatialCoords(sfe), ignore_attr = "dimnames")
})
