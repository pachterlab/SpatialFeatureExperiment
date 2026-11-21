library(SpatialExperiment)
library(SingleCellExperiment)
library(DropletUtils)
library(sf)
library(SFEData)
library(terra)
library(VisiumIO)

if (!dir.exists("ob")) dir.create("ob", recursive = TRUE)
mat_fn <- file.path("ob", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
                  destfile = file.path("ob", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("ob", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz",
                  destfile = file.path("ob", "spatial.tar.gz"))
    untar(file.path("ob", "spatial.tar.gz"), exdir = "ob")
}

sce <- read10xCounts(file.path("ob", "filtered_feature_bc_matrix.h5"))
colnames(sce) <- sce$Barcode

# Read in spatial coordinates
coords <- read.csv(file.path("ob", "spatial", "tissue_positions.csv"), row.names = 1)
obs <- intersect(colnames(sce), rownames(coords))
sce <- sce[,obs]
coords <- coords[obs,]
int_colData(sce)$spatialCoords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
spe <- as(sce, "SpatialExperiment")

test_that("Convert SPE and SCE to SFE, no images", {
    sfe1 <- as(spe, "SpatialFeatureExperiment")
    sfe2 <- as(sce, "SpatialFeatureExperiment")
    sfe3 <- toSpatialFeatureExperiment(spe)
    #colData(sce) <- cbind(colData(sce), coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
    sfe4 <- toSpatialFeatureExperiment(sce, spatialCoordsNames = NULL)
    sfe <- read10xVisiumSFE(dirs = "ob", data = "filtered", flip = "image", add_centroids = TRUE)
    centroids_check <- st_geometry(centroids(sfe))

    expect_s4_class(sfe1, "SpatialFeatureExperiment")
    expect_equal(sfe2, sfe1)
    expect_equal(sfe3, sfe1)
    expect_equal(sfe4, sfe1)
    expect_equal(unname(unlist(st_equals(st_geometry(centroids(sfe1)), centroids_check))),
                 seq_len(ncol(sfe1)))
    expect_equal(spatialCoords(sfe1), spatialCoords(sfe), ignore_attr = "dimnames")
})

test_that("Convert SPE to SFE, loaded images", {
    spe <- TENxVisium(spacerangerOut = "ob", processing = "filtered", format = "h5",
                      images = c("lowres"), loadImage = TRUE) |> import()
    sfe <- toSpatialFeatureExperiment(spe)
    img1 <- getImg(spe)
    img2 <- getImg(sfe)
    expect_s4_class(img2, "SpatRasterImage")
    expect_equal(dim(img1), dim(img2)[1:2])
    v1 <- col2rgb(imgRaster(img1))
    v2 <- terra::values(img2)
    v2 <- t(v2)
    dimnames(v1) <- dimnames(v2) <- NULL
    expect_equal(v1, v2)
    
    # Image is properly aligned in the SFE
    v <- terra::extract(terra::mean(img2), centroids(sfe) |> vect())
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mean) < -0.4)
})

test_that("Convert SPE to SFE, stored images", {
    spe <- TENxVisium(spacerangerOut = "ob", processing = "filtered", format = "h5",
                      images = c("lowres")) |> import()
    sfe <- toSpatialFeatureExperiment(spe)
    img1 <- getImg(spe)
    img2 <- getImg(sfe)
    expect_s4_class(img2, "SpatRasterImage")
    expect_equal(dim(img1), dim(img2)[1:2])
    v1 <- col2rgb(imgRaster(img1))
    v2 <- terra::values(img2)
    v2 <- t(v2)
    dimnames(v1) <- dimnames(v2) <- NULL
    expect_equal(v1, v2)

    # Image is properly aligned in the SFE
    v <- terra::extract(terra::mean(img2), centroids(sfe) |> vect())
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mean) < -0.4)
})

test_that("Convert SPE to SFE, with SpatRaster image", {
    spe <- TENxVisium(spacerangerOut = "ob", processing = "filtered", format = "h5",
                      images = c("lowres")) |> import()
    suppressWarnings(Img(spe, image_id = "lowres") <- rast(imgSource(getImg(spe))) |> SpatRasterImage())
    sfe <- toSpatialFeatureExperiment(spe)
    expect_true(is(getImg(sfe), "SpatRaster"))
    expect_equal(getImg(spe), getImg(sfe))
})

test_that("Convert SPE to SFE, when the SPE has BioFormatsImage", {
    fp <- tempfile()
    xenium <- XeniumOutput("v1")
    sfe <- readXenium(xenium)
    spe <- as(sfe, "SpatialExperiment")
    sfe2 <- toSpatialFeatureExperiment(spe)
    expect_s4_class(getImg(sfe2), "BioFormatsImage")
})
