outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
bc_flou1 <- read.csv(file.path(outdir, "sample01", "outs", "spatial",
                               "barcode_fluorescence_intensity.csv"))
sp_enr1 <- read.csv(file.path(outdir, "sample01", "outs", "spatial",
                              "spatial_enrichment.csv"))
pos1 <- read.csv(file.path(outdir, "sample01", "outs", "spatial",
                           "tissue_positions.csv"))

bc_flou2 <- read.csv(file.path(outdir, "sample02", "outs", "spatial",
                               "barcode_fluorescence_intensity.csv"))
sp_enr2 <- read.csv(file.path(outdir, "sample02", "outs", "spatial",
                              "spatial_enrichment.csv"))
pos2 <- read.csv(file.path(outdir, "sample02", "outs", "spatial",
                           "tissue_positions.csv"))

samples <- file.path(outdir, paste0("sample0", 1:2))

rd1 <- sp_enr1[,4:9]
rd2 <- sp_enr2[,4:9]
names(rd1) <- paste(names(rd1), "sample01", sep = "_")
names(rd2) <- paste(names(rd2), "sample02", sep = "_")

rd_expect <- cbind(Feature.Type = sp_enr1$Feature.Type, rd1, rd2)

cd_expect <- rbind(bc_flou1, bc_flou2)[, 3:8]

test_that("Correctly read Space Ranger output", {
    sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
    # Very basic one
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(sampleIDs(sfe), c("sample01", "sample02"))
    expect_equal(colGeometryNames(sfe), "spotPoly")
    expect_equal(colGraphNames(sfe, "sample01"), "visium")
    expect_equal(colGraphNames(sfe, "sample02"), "visium")
    expect_equal(as.data.frame(rowData(sfe)[,-1]), rd_expect,
                 ignore_attr = "row.names")
    expect_equal(as.data.frame(colData(sfe)[,5:10]), cd_expect,
                 ignore_attr = "row.names")
})

test_that("Correctly add visium graph when there's 1 sample", {
    sfe <- read10xVisiumSFE(samples[1], type = "sparse", data = "filtered")
    expect_equal(colGraphNames(sfe, "sample01"), "visium")
})

# Need uncropped image
if (!dir.exists("ob")) dir.create(file.path("ob", "outs"), recursive = TRUE)
mat_fn <- file.path("ob", "outs", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
                  destfile = file.path("ob", "outs", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("ob", "outs", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz",
                  destfile = file.path("ob", "outs", "spatial.tar.gz"))
    untar(file.path("ob", "outs", "spatial.tar.gz"), exdir = file.path("ob", "outs"))
}

if (!dir.exists("kidney")) dir.create(file.path("kidney", "outs"), recursive = TRUE)
mat_fn <- file.path("kidney", "outs", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_filtered_feature_bc_matrix.h5",
                  destfile = file.path("kidney", "outs", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("kidney", "outs", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_spatial.tar.gz",
                  destfile = file.path("kidney", "outs", "spatial.tar.gz"))
    untar(file.path("kidney", "outs", "spatial.tar.gz"), exdir = file.path("kidney", "outs"))
}

library(terra)
library(sf)
library(SingleCellExperiment)
library(SpatialExperiment)

test_that("Image is properly aligned in pixel space", {
    sfe <- read10xVisiumSFE("ob")
    expect_equal(unit(sfe), "full_res_image_pixel")
    cg <- spotPoly(sfe)
    cg$nCounts <- Matrix::colSums(counts(sfe))
    cg$geometry <- st_centroid(cg$geometry)
    img_lo <- getImg(sfe, image_id = "lowres") |> imgRaster()
    img_lo <- terra::mean(img_lo)
    v_lo <- terra::extract(img_lo, cg)
    # This test only works for this tissue for filtered data
    expect_true(abs(cor(cg$nCounts, v_lo$mean)) > 0.4)
    img_hi <- getImg(sfe, image_id = "hires") |> imgRaster()
    img_hi <- terra::mean(img_hi)
    v_hi <- terra::extract(img_hi, cg)
    expect_true(abs(cor(cg$nCounts, v_hi$mean)) > 0.4)
    bbox_cg <- st_as_sfc(st_bbox(cg))
    bbox_img_lo <- st_as_sfc(st_bbox(as.vector(ext(img_lo))))
    bbox_img_hi <- st_as_sfc(st_bbox(as.vector(ext(img_hi))))
    expect_equal(bbox_img_lo, bbox_img_hi)
    expect_true(st_covered_by(bbox_cg, bbox_img_lo, sparse = FALSE))
    expect_true(st_area(bbox_cg)/st_area(bbox_img_lo) > 0.1)
})

test_that("Read when one out of multiple images are desired", {
    sfe <- read10xVisiumSFE("ob", images = "lowres")
    expect_equal(nrow(imgData(sfe)), 1L)
    expect_equal(imgData(sfe)$image_id, "lowres")
})

test_that("Image is properly aligned in micron space", {
    sfe <- read10xVisiumSFE("ob", unit = "micron")
    expect_equal(unit(sfe), "micron")
    cg <- spotPoly(sfe)
    cg$nCounts <- Matrix::colSums(counts(sfe))
    cg$geometry <- st_centroid(cg$geometry)
    img_lo <- getImg(sfe, image_id = "lowres") |> imgRaster()
    img_lo <- terra::mean(img_lo)
    v_lo <- terra::extract(img_lo, cg)
    expect_true(abs(cor(cg$nCounts, v_lo$mean)) > 0.4)
    img_hi <- getImg(sfe, image_id = "hires") |> imgRaster()
    img_hi <- terra::mean(img_hi)
    v_hi <- terra::extract(img_hi, cg)
    expect_true(abs(cor(cg$nCounts, v_hi$mean)) > 0.4)
    bbox_cg <- st_as_sfc(st_bbox(cg))
    bbox_img_lo <- st_as_sfc(st_bbox(as.vector(ext(img_lo))))
    bbox_img_hi <- st_as_sfc(st_bbox(as.vector(ext(img_hi))))
    expect_equal(bbox_img_lo, bbox_img_hi)
    expect_true(st_covered_by(bbox_cg, bbox_img_lo, sparse = FALSE))
    expect_true(st_area(bbox_cg)/st_area(bbox_img_lo) > 0.1)
    expect_true(all(st_coordinates(bbox_img_lo) < 1e4))
    expect_true(all(st_coordinates(bbox_cg) < 1e4))
})

test_that("Micron spot spacing works when there're singletons", {
    sfe <- read10xVisiumSFE("kidney", unit = "micron", zero.policy = TRUE)
    expect_equal(unit(sfe), "micron")
})

dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")

if (!dir.exists("vizgen")) {
    file.copy(dir_use, ".", recursive = TRUE)
}
dir_use <- "vizgen"
test_that("readVizgen flip geometry, use cellpose", {
    sfe <- readVizgen(dir_use, z = 0L, use_cellpose = TRUE, image = "PolyT",
                      flip = "geometry", min_area = 15)
    expect_equal(unit(sfe), "micron")
    expect_equal(imgData(sfe)$image_id, "PolyT_z0")
    img <- imgRaster(getImg(sfe))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mosaic_PolyT_z0) > 0.4)
    # Make sure both segmentations and centroids are flipped
    hulls <- st_convex_hull(cellSeg(sfe))
    expect_true(all(vapply(seq_len(nrow(cg)), function(i) {
        st_covered_by(cg[i,], hulls[i,], sparse = FALSE)[1,1]
    }, FUN.VALUE = logical(1))))
    # Make sure that cells that are too small are removed
    cg <- cellSeg(sfe)
    areas <- st_area(cg)
    expect_true(all(areas > 15))
})

test_that("readVizgen flip geometry, don't use cellpose", {
    sfe <- readVizgen(dir_use, z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "geometry")
    expect_equal(unit(sfe), "micron")
    img <- imgRaster(getImg(sfe))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mosaic_PolyT_z0) > 0.4)
    hulls <- st_convex_hull(cellSeg(sfe))
    expect_true(all(vapply(seq_len(nrow(cg)), function(i) {
        st_covered_by(cg[i,], hulls[i,], sparse = FALSE)[1,1]
    }, FUN.VALUE = logical(1))))
})

test_that("readVizgen flip image", {
    sfe <- readVizgen(dir_use, z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image")
    expect_equal(unit(sfe), "micron")
    img <- imgRaster(getImg(sfe))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mosaic_PolyT_z0) > 0.4)
})

test_that("readVizgen don't flip image when image is too large", {
    expect_error(readVizgen(dir_use, z = 0L, use_cellpose = FALSE, image = "PolyT",
                            flip = "image", max_flip = "0.02 TB"),
                 "max_flip must be in either MB or GB")
    sfe <- readVizgen(dir_use, z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image", max_flip = "0.02 MB")
    suppressWarnings(img_orig <- rast(file.path(dir_use, "images", "mosaic_PolyT_z0.tif")))
    img <- imgRaster(getImg(sfe))
    # Make sure image was not flipped
    expect_equal(terra::values(img), terra::values(img_orig))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mosaic_PolyT_z0) > 0.4)
})

test_that("Don't flip image if it's GeoTIFF", {
    sfe <- readVizgen(dir_use, z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image")
    terra::writeRaster(imgRaster(getImg(sfe)),
                       filename = file.path("vizgen", "images", "mosaic_DAPI_z0.tif"),
                       overwrite = TRUE)
    sfe2 <- readVizgen(dir_use, z = 0L, use_cellpose = FALSE, image = "DAPI",
                       flip = "image")
    expect_equal(terra::values(imgRaster(getImg(sfe))), terra::values(imgRaster(getImg(sfe2))))
    file.remove(file.path("vizgen", "images", "mosaic_DAPI_z0.tif"))
})

test_that("Errors and warnings", {
    expect_warning(sfe <- readVizgen(dir_use, z = 0L, image = "DAPI", use_cellpose = FALSE),
                   "don't exist")
    expect_equal(nrow(imgData(sfe)), 0L)
    expect_error(readVizgen(dir_use, z = 7L, image = "PolyT", use_cellpose = FALSE),
                 "z must be beween 0 and 6")
})

# Make toy examples of multiple pieces
parq <- sfarrow::st_read_parquet(file.path(dir_use, "cellpose_micron_space.parquet"))
parq2 <- parq[1:4,]
# One large piece and one small piece
large <- list(matrix(c(2500, 0,
                       2510, 0,
                       2510, 10,
                       2500, 10,
                       2500, 0), ncol = 2, byrow = TRUE))
small <- list(matrix(c(2515, 0,
                       2516, 0,
                       2516, 1,
                       2515, 0), ncol = 2, byrow = TRUE))
small2 <- list(small[[1]] + 5)
large2 <- list(large[[1]] * 0.9 + 20)
large_small <- st_multipolygon(list(large, small))
large_g <- st_multipolygon(list(large))
small_small <- st_multipolygon(list(small, small2))
large_large <- st_multipolygon(list(large, large2))

new_geo <- st_sfc(large_g, large_small, small_small, large_large)
parq2$Geometry <- new_geo
if (!dir.exists("multi")) {
    dir.create("multi")
    file.copy("vizgen", "multi", recursive = TRUE)
}

dir_use <- file.path("multi", "vizgen")
file.remove(file.path(dir_use, "cellpose_micron_space.parquet"))
suppressWarnings(sfarrow::st_write_parquet(parq2, file.path(dir_use, "cellpose_micron_space.parquet")))
test_that("Deal with multiple pieces, remove pieces that are too small", {
    w <- capture_warnings(sfe <- readVizgen(dir_use, z = 0L, image = "PolyT"))
    expect_match(w, "Sanity check", all = FALSE)
    expect_match(w, "The largest piece is kept", all = FALSE)
    cg <- cellSeg(sfe)
    expect_equal(st_geometry_type(cg, by_geometry = "FALSE") |> as.character(), "POLYGON")
    expect_equal(colnames(sfe), parq2$EntityID[c(1,2,4)])
    areas <- st_area(cg)
    expect_true(all(vapply(areas, all.equal, target = st_area(large_g),
                           FUN.VALUE = logical(1))))
})

# Read all z-planes
# When some cells are removed because they're too small
# Message when removing empty geometries when reading Vizgen
test_that("Read multiple z-planes for Vizgen", {
# Need to make toy data for multiple z-planes
})

test_that("Read different versions of Vizgen data", {
    # Version with HDF5 already tested

})

test_that("Read MERFISH transcript spots into rowGeometries", {

})

test_that("Read MERFISH transcript spots into colGeometries", {

})

test_that("Read MERFISH transcript spots into imgData", {

})
