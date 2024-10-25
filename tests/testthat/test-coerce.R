library(SpatialExperiment)
library(SingleCellExperiment)
library(DropletUtils)
library(sf)

dir <- system.file("extdata/sample01", package = "SpatialFeatureExperiment")

sce <- read10xCounts(file.path(dir, "outs", "filtered_feature_bc_matrix"))

# Read in spatial coordinates
coords <- read.csv(file.path(dir, "outs", "spatial", "tissue_positions.csv"))
int_colData(sce)$spatialCoords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
spe <- as(sce, "SpatialExperiment")

test_that("Convert SPE and SCE to SFE, no images", {
    sfe1 <- as(spe, "SpatialFeatureExperiment")
    sfe2 <- as(sce, "SpatialFeatureExperiment")
    sfe3 <- toSpatialFeatureExperiment(spe)
    colData(sce) <- cbind(colData(sce), coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
    sfe4 <- toSpatialFeatureExperiment(sce, spatialCoordsNames = c("pxl_col_in_fullres",
                                                                   "pxl_row_in_fullres"))
    sfe <- read10xVisiumSFE(dir, type = "sparse", data = "filtered", flip = "image")
    centroids_check <- st_centroid(st_geometry(spotPoly(sfe)))

    expect_s4_class(sfe1, "SpatialFeatureExperiment")
    expect_equal(sfe2, sfe1)
    expect_equal(sfe3, sfe1)
    expect_equal(sfe4, sfe1)
    expect_equal(st_geometry(centroids(sfe1)), centroids_check, ignore_attr = "row.names")
    expect_equal(spatialCoords(sfe1), spatialCoords(sfe), ignore_attr = "dimnames")
})

test_that("Convert SPE to SFE, loaded images", {
    spe <- read10xVisium(dir, load = TRUE, type = "sparse")
    sfe <- toSpatialFeatureExperiment(spe)
    img1 <- getImg(spe)
    img2 <- getImg(sfe)
    expect_s4_class(img2, "SpatRasterImage")
    expect_equal(dim(img1), dim(img2)[1:2])
    v1 <- col2rgb(imgRaster(img1))
    v2 <- terra::values(imgRaster(img2))
    v2 <- t(v2)
    dimnames(v1) <- dimnames(v2) <- NULL
    expect_equal(v1, v2)

    bbox <- st_bbox(centroids(sfe))
    bbox_img <- as.vector(ext(imgRaster(img2))) # That the image is properly scaled
    diffs1 <- bbox[3:4] - bbox[1:2]
    diffs2 <- bbox_img[c(2,4)] - bbox_img[c(1,3)]
    expect_true(all(diffs1 / diffs2 > (1-1/min(dim(img1)))))
})

test_that("Convert SPE to SFE, stored images", {
    spe <- read10xVisium(dir, load = FALSE, type = "sparse")
    sfe <- toSpatialFeatureExperiment(spe)
    img1 <- getImg(spe)
    img2 <- getImg(sfe)
    expect_s4_class(img2, "SpatRasterImage")
    expect_equal(dim(imgRaster(img1)), dim(img2)[1:2])
    v1 <- col2rgb(imgRaster(img1))
    v2 <- terra::values(imgRaster(img2))
    v2 <- t(v2)
    dimnames(v1) <- dimnames(v2) <- NULL
    expect_equal(v1, v2)

    bbox <- st_bbox(centroids(sfe))
    bbox_img <- as.vector(ext(imgRaster(img2)))
    diffs1 <- bbox[3:4] - bbox[1:2]
    diffs2 <- bbox_img[c(2,4)] - bbox_img[c(1,3)]
    expect_true(all(diffs1 / diffs2 > (1-1/min(dim(img1)))))
})
