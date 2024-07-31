library(SFEData)
library(sfarrow)
library(S4Vectors)
# Read Visium=============
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
# Read Vizgen MERFISH==============
test_that("readVizgen flip geometry, use cellpose", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    expect_message(sfe <- readVizgen(dir_use, z = 3L, use_cellpose = TRUE,
                                     flip = "geometry", min_area = 15),
                   "with area less than 15")
    expect_equal(unit(sfe), "micron")
    expect_equal(imgData(sfe)$image_id,
                 paste0(c(paste0("Cellbound", 1:3), "DAPI", "PolyT"),
                       "_z3"))
    img <- imgRaster(getImg(sfe, image_id = "PolyT_z3"))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    # Shouldn't be many cells in that empty space if properly aligned
    expect_true(sum(v$mosaic_PolyT_z3 < 30, na.rm = TRUE) < 10)
    # Make sure both segmentations and centroids are flipped
    hulls <- st_convex_hull(cellSeg(sfe))
    # Geometries were literally cropped so centroids may be outside cropped cells
    bbox_use <- st_bbox(hulls) |> st_as_sfc() |> st_cast("LINESTRING")
    boundary <- st_touches(hulls, bbox_use)
    inds <- which(lengths(boundary) == 0L)
    hulls2 <- hulls[inds,]
    cg2 <- cg[inds,]
    expect_true(all(vapply(seq_len(nrow(cg2)), function(i) {
        st_covered_by(cg2[i,], hulls2[i,], sparse = FALSE)[1,1]
    }, FUN.VALUE = logical(1))))
    # Make sure that cells that are too small are removed
    cg <- cellSeg(sfe)
    areas <- st_area(cg)
    expect_true(all(areas > 15))
    unlink(dir_use, recursive = TRUE)
})

test_that("readVizgen write parquet file after reading hdf5", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    file.remove(file.path(dir_use, "cell_boundaries.parquet"))
    sfe <- readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "PolyT")
    expect_true(file.exists(file.path(dir_use,"hdf5s_micron_space.parquet")))
    # Second time reading
    expect_message(readVizgen(dir_use, z = 3L, image = "PolyT"),
                   "processed hdf5 files will be used")
    unlink(dir_use, recursive = TRUE)
})

test_that("readVizgen flip geometry, don't use cellpose", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    sfe <- readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "PolyT",
                      flip = "geometry")
    expect_equal(unit(sfe), "micron")
    img <- imgRaster(getImg(sfe))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    expect_true(sum(v$mosaic_PolyT_z3 < 30, na.rm = TRUE) < 10)
    # Make sure both segmentations and centroids are flipped
    hulls <- st_convex_hull(cellSeg(sfe))
    # Geometries were literally cropped so centroids may be outside cropped cells
    bbox_use <- st_bbox(hulls) |> st_as_sfc() |> st_cast("LINESTRING")
    boundary <- st_touches(hulls, bbox_use)
    inds <- which(lengths(boundary) == 0L)
    hulls2 <- hulls[inds,]
    cg2 <- cg[inds,]
    expect_true(all(vapply(seq_len(nrow(cg2)), function(i) {
        st_covered_by(cg2[i,], hulls2[i,], sparse = FALSE)[1,1]
    }, FUN.VALUE = logical(1))))
    unlink(dir_use, recursive = TRUE)
})

test_that("readVizgen flip image", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    sfe <- readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image")
    expect_equal(unit(sfe), "micron")
    img <- imgRaster(getImg(sfe))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    expect_true(sum(v$mosaic_PolyT_z3 < 30, na.rm = TRUE) < 10)
    unlink(dir_use, recursive = TRUE)
})

test_that("readVizgen don't flip image when image is too large", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    expect_error(readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "PolyT",
                            flip = "image", max_flip = "0.02 TB"),
                 "max_flip must be in either MB or GB")
    sfe <- readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image", max_flip = "0.02 MB")
    suppressWarnings(img_orig <- rast(file.path(dir_use, "images", "mosaic_PolyT_z3.tif")))
    img <- imgRaster(getImg(sfe))
    # Make sure image was not flipped
    expect_equal(terra::values(img), terra::values(img_orig))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    expect_true(sum(v$mosaic_PolyT_z3 < 30, na.rm = TRUE) < 10)
    unlink(dir_use, recursive = TRUE)
})

test_that("Don't flip image if it's GeoTIFF", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    sfe <- readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image")
    terra::writeRaster(imgRaster(getImg(sfe)),
                       filename = file.path(dir_use, "images", "mosaic_DAPI_z3.tif"),
                       overwrite = TRUE)
    sfe2 <- readVizgen(dir_use, z = 3L, use_cellpose = FALSE, image = "DAPI",
                       flip = "image")
    expect_equal(terra::values(imgRaster(getImg(sfe))), terra::values(imgRaster(getImg(sfe2))))
    unlink(dir_use, recursive = TRUE)
})

test_that("Errors and warnings in readVizgen", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    expect_warning(sfe <- readVizgen(dir_use, z = 1L, image = "PolyT"),
                   "don't exist")
    expect_equal(nrow(imgData(sfe)), 0L)
    expect_error(readVizgen(dir_use, z = 7L, image = "PolyT"),
                 "z must be beween 0 and 6")
    unlink(dir_use, recursive = TRUE)
})

# Make toy examples of multiple pieces
fp <- tempdir()
dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
parq <- st_read_parquet(file.path(dir_use, "cell_boundaries.parquet"))
unlink(dir_use, recursive = TRUE)
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

test_that("Deal with multiple pieces, remove pieces that are too small", {
    fp <- tempdir()
    parq2 <- parq[1:4,]
    new_geo <- st_sfc(large_g, large_small, small_small, large_large)
    parq2$Geometry <- new_geo
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "multi"))
    file.remove(file.path(dir_use, "cell_boundaries.parquet"))
    suppressWarnings(st_write_parquet(parq2, file.path(dir_use, "cell_boundaries.parquet")))

    w <- capture_warnings(sfe <- readVizgen(dir_use, z = 3L, image = "PolyT"))
    expect_match(w, "Sanity check", all = FALSE)
    expect_match(w, "The largest piece is kept", all = FALSE)
    cg <- cellSeg(sfe)
    expect_equal(st_geometry_type(cg, by_geometry = "FALSE") |> as.character(), "POLYGON")
    expect_equal(colnames(sfe), as.character(parq2$EntityID[c(1,2,4)]))
    areas <- st_area(cg)
    expect_true(all(vapply(areas, all.equal, target = st_area(large_g),
                           FUN.VALUE = logical(1))))
    unlink(dir_use, recursive = TRUE)
})

test_that("No polygons left", {
    # Like they're all too small, or when the polygon file is empty, unlikely
    # but otherwise we get a mysterious error
    fp <- tempdir()
    parq2 <- parq[1:2,]
    small <- st_polygon(small)
    new_geo <- st_sfc(small, small)
    parq2$Geometry <- new_geo
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "small"))

    file.remove(file.path(dir_use, "cell_boundaries.parquet"))
    suppressWarnings(st_write_parquet(parq2, file.path(dir_use, "cell_boundaries.parquet")))

    expect_error(readVizgen(dir_use, z = 3L, image = "PolyT"),
                 "No polygons left after filtering")
    unlink(dir_use, recursive = TRUE)
})

test_that("Image can be named polyT in older version", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "image"))

    file.rename(file.path(dir_use, "images", "mosaic_PolyT_z3.tif"),
                file.path(dir_use, "images", "mosaic_polyT_z3.tif"))
    expect_no_warning(readVizgen(dir_use, z = 3L, image = "PolyT"))
    unlink(dir_use, recursive = TRUE)
})

test_that("Version with Cellpose directory", {
    fp <- tempdir()
    dir_use <- VizgenOutput("cellpose", file_path = file.path(fp, "cellpose"))
    # Cropped geometry has 2 pieces in one cell
    suppressWarnings(sfe <- readVizgen(dir_use, z = 3L, use_cellpose = TRUE,
                                       flip = "geometry", min_area = 15))
    expect_equal(unit(sfe), "micron")
    expect_equal(imgData(sfe)$image_id,
                 paste0(c(paste0("Cellbound", 1:3), "DAPI", "PolyT"),
                        "_z3"))
    img <- imgRaster(getImg(sfe, image_id = "PolyT_z3"))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    # Shouldn't be many cells in that empty space if properly aligned
    expect_true(sum(v$mosaic_PolyT_z3 < 30, na.rm = TRUE) < 10)
    # Make sure both segmentations and centroids are flipped
    hulls <- st_convex_hull(cellSeg(sfe))
    # Geometries were literally cropped so centroids may be outside cropped cells
    bbox_use <- st_bbox(hulls) |> st_as_sfc() |> st_cast("LINESTRING")
    boundary <- st_touches(hulls, bbox_use)
    inds <- which(lengths(boundary) == 0L)
    hulls2 <- hulls[inds,]
    cg2 <- cg[inds,]
    expect_true(all(vapply(seq_len(nrow(cg2)), function(i) {
        st_covered_by(cg2[i,], hulls2[i,], sparse = FALSE)[1,1]
    }, FUN.VALUE = logical(1))))
    # Make sure that cells that are too small are removed
    cg <- cellSeg(sfe)
    areas <- st_area(cg)
    expect_true(all(areas > 15))
    unlink(dir_use, recursive = TRUE)
})

test_that("Message when removing empty polygons", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "empty"))

    parq <- st_read_parquet(file.path(dir_use, "cell_boundaries.parquet"))
    parq$Geometry[[1]] <- st_polygon()
    file.remove(file.path(dir_use, "cell_boundaries.parquet"))
    suppressWarnings(st_write_parquet(parq, file.path(dir_use, "cell_boundaries.parquet")))

    expect_message(sfe <- readVizgen(dir_use, z = 3L, image = "PolyT"),
                   "..removing 1 empty polygons")
    unlink(dir_use, recursive = TRUE)
})

test_that("Read all z-planes for Vizgen", {
    fp <- tempdir()
    # Only affecting images
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "all"))
    sfe <- readVizgen(dir_use, z = "all", image = "DAPI")
    expect_equal(imgData(sfe)$image_id, paste0("DAPI_z", 1:3))
    unlink(dir_use, recursive = TRUE)
})

test_that("Error message when multiple parquet files are found in readVizgen", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    file.copy(file.path(dir_use, "cell_boundaries.parquet"),
              file.path(dir_use, "cellpose_micron_space.parquet"))
    m <- capture_messages(sfe <- readVizgen(dir_use, z = "all", image = "DAPI"))
    expect_match(m, " `.parquet` files exist:", all = FALSE)
    expect_match(m, "using ->", all = FALSE)

    file.rename(file.path(dir_use, "cellpose_micron_space.parquet"),
                file.path(dir_use, "cool_cell_boundaries.parquet"))
    expect_error(sfe <- readVizgen(dir_use, z = "all", image = "DAPI"),
                 "only 1 `.parquet` file can be read")
    unlink(dir_use, recursive = TRUE)
})

test_that("Make cell bounding boxes when segmentation is absent", {
    fp2 <- tempdir()
    dir_use <- VizgenOutput("cellpose", file_path = file.path(fp, "vizgen_test"))
    file.remove(file.path(dir_use, "Cellpose", "cellpose_micron_space.parquet"))
    expect_message(
        expect_warning(sfe <- readVizgen(dir_use, use_bboxes = TRUE),
                       "No '.hdf5' files present"),
        "Creating bounding boxes")
    expect_equal(colGeometryNames(sfe), c("centroids", "cell_bboxes"))
    unlink(dir_use, recursive = TRUE)
})

# Read CosMX===================
test_that("readCosMX, not reading transcript spots", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))
    sfe <- readCosMX(dir_use, z = 1L)
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), c("centroids", "cellSeg"))
    expect_equal(st_geometry_type(cellSeg(sfe), by_geometry = FALSE) |> as.character(),
                 "POLYGON")
    expect_equal(dim(sfe), c(960, 27))
    # parquet file for cell polygons written
    expect_true(file.exists(file.path(dir_use, "cell_boundaries_sf.parquet")))
    time_note <- Sys.time()
    # Second time reading
    sfe <- readCosMX(dir_use, z = 1L)
    time_file <- file.info(file.path(dir_use, "cell_boundaries_sf.parquet"))$ctime
    expect_true(time_file < time_note)
    unlink(dir_use, recursive = TRUE)
})

test_that("readCosMX, reading spots, 1 z-plane", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))
    sfe <- readCosMX(dir_use, z = 1L, add_molecules = TRUE)
    fn <- file.path(dir_use, "tx_spots.parquet")
    expect_true(file.exists(fn))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(st_geometry_type(txSpots(sfe), FALSE) |> as.character(),
                 "MULTIPOINT")
    expect_null(st_z_range(txSpots(sfe)))

    # Reloading the second time reading
    time_note <- Sys.time() # Check the files weren't written again
    sfe <- readCosMX(dir_use, z = 1L, add_molecules = TRUE)
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(st_geometry_type(txSpots(sfe), FALSE) |> as.character(),
                 "MULTIPOINT")
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    unlink(dir_use, recursive = TRUE)
})

test_that("readCosMX, 2 z-planes, split z, not splitting by cell compartments", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))

    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     z_option = "split")
    d <- file.path(dir_use, "tx_spots")
    expect_true(dir.exists(d))
    fns_expect <- paste0("tx_spots_z", 0:1, ".parquet")
    expect_equal(list.files(d), fns_expect)
    expect_equal(rowGeometryNames(sfe), paste0("tx_spots_z", 0:1))
    expect_null(st_z_range(rowGeometry(sfe, "tx_spots_z0")))

    # Reloading the second time reading
    # Both z-planes
    time_note <- Sys.time()
    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     z_option = "split")
    expect_null(st_z_range(rowGeometry(sfe, "tx_spots_z0")))
    fn <- file.path(d, "tx_spots_z0.parquet")
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    # Only read one of the z-planes
    sfe <- readCosMX(dir_use, z = 0L, add_molecules = TRUE)
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    unlink(dir_use, recursive = TRUE)
})

test_that("readCosMX, don't split z, don't split cell compartments", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))

    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     z_option = "3d")
    fn <- file.path(dir_use, "tx_spots.parquet")
    expect_true(file.exists(fn))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(unclass(st_z_range(txSpots(sfe))), c(zmin = 0, zmax = 1),
                 ignore_attr = "crs")

    unlink(dir_use, recursive = TRUE)
})

test_that("readCosMX, split z, split cell compartments", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))

    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     split_cell_comps = TRUE, z_option = "split")
    d <- file.path(dir_use, "tx_spots")
    expect_true(dir.exists(d))
    comps <- c("Nuclear", "None", "Membrane", "Cytoplasm")
    combs <- expand.grid(compartment = comps, z = 0:1, stringsAsFactors = FALSE)
    rgns <- paste0(combs$compartment, "_z", combs$z)
    fns <- paste0(rgns, ".parquet")

    expect_setequal(rowGeometryNames(sfe), rgns)
    expect_setequal(list.files(d), fns)

    # Reloading the second time reading
    time_note <- Sys.time()
    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     split_cell_comps = TRUE, z_option = "split")
    time_check <- file.info(file.path(d, fns[1]))$mtime
    expect_true(time_note > time_check)
    expect_setequal(rowGeometryNames(sfe), rgns)

    # Only read one of the z-planes
    sfe <- readCosMX(dir_use, z = 0L, add_molecules = TRUE,
                     split_cell_comps = TRUE)
    time_check <- file.info(file.path(d, fns[1]))$mtime
    expect_true(time_note > time_check)
    expect_setequal(rowGeometryNames(sfe), rgns[1:4])

    unlink(dir_use, recursive = TRUE)
})

# Read Xenium XOA v1================
# Flip image
test_that("readXenium, XOA v1", {
    library(RBioFormats)
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    # Weirdly the first time I get the null pointer error
    try(m <- read.omexml(file.path(fn, "morphology_focus.ome.tif")))
    sfe <- readXenium(fn, add_molecules = TRUE)
    # Basic stuff
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), c("centroids", "cellSeg", "nucSeg"))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(as.character(st_geometry_type(SpatialFeatureExperiment::centroids(sfe), FALSE)), "POINT")
    expect_equal(as.character(st_geometry_type(cellSeg(sfe), FALSE)), "POLYGON")
    expect_equal(as.character(st_geometry_type(txSpots(sfe), FALSE)), "MULTIPOINT")
    expect_equal(imageIDs(sfe), c("morphology_focus", "morphology_mip"))
    expect_s4_class(getImg(sfe, image_id = "morphology_focus"), "BioFormatsImage")
    expect_s4_class(getImg(sfe, image_id = "morphology_mip"), "BioFormatsImage")
    expect_equal(dim(getImg(sfe, image_id = "morphology_focus"))[["C"]], 1L)

    # That things are aligned
    bbox_rg <- st_bbox(txSpots(sfe)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))

    img <- toExtImage(getImg(sfe), resolution = 1L)
    mask <- img > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe)))
    expect_true(mean(v$lyr.1) > 0.9)
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v1, image not found", {
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    file.remove(file.path(fn, "morphology_focus.ome.tif"))
    expect_warning(sfe <- readXenium(fn, add_molecules = TRUE),
                   "or have non-standard file name")
    expect_equal(imageIDs(sfe), "morphology_mip")
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v1, use parquet files, with annoying arrow raw bytes", {
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    file.remove(file.path(fn, "cell_boundaries.csv.gz"))
    file.remove(file.path(fn, "nucleus_boundaries.csv.gz"))
    file.rename(file.path(fn, "cell_boundaries_binary.parquet"),
                file.path(fn, "cell_boundaries.parquet"))
    file.rename(file.path(fn, "nucleus_boundaries_binary.parquet"),
                file.path(fn, "nucleus_boundaries.parquet"))
    file.remove(list.files(fn, pattern = "nobinary.parquet", full.names = TRUE))
    expect_message(sfe <- readXenium(fn), "Converting columns with raw bytes")
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v1 when only cell but not nuclei segmentation is available", {
    # Since the `polys` object should be a list even if there's only one of cell or nuclei
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    file.remove(list.files(fn, "cell_boundaries\\.*", full.names = TRUE))
    expect_warning(sfe <- readXenium(fn), 'No `cell_boundaries` file is available')
    expect_equal(colGeometryNames(sfe), c("centroids", "nucSeg"))
    expect_warning(sfe2 <- readXenium(fn, segmentations = "cell"),
                    'No `cell_boundaries` file is available')
    expect_equal(colGeometryNames(sfe2), "centroids")
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v1 read the output _sf.parquet next time", {
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    sfe <- readXenium(fn)
    expect_true(file.exists(file.path(fn, "cell_boundaries_sf.parquet")))
    expect_true(file.exists(file.path(fn, "nucleus_boundaries_sf.parquet")))
    time_note <- Sys.time()
    sfe <- readXenium(fn)
    time_check <- file.info(file.path(fn, "cell_boundaries_sf.parquet"))$mtime
    expect_true(time_note > time_check)
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v1 read cell metadata parquet when csv is absent", {
    skip_on_bioc() # zstd error, only on mac, don't know what to do
    # Not sure if it's about Xenium v1 or the way the subset was written
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    file.remove(file.path(fn, "cells.csv.gz"))
    expect_message(sfe <- readXenium(fn), ">>> Reading cell metadata -> `cells.parquet`")
    unlink(fn, recursive = TRUE)
}) # Would be nice though to use csv if arrow is not installed in case the user
# has trouble installing arrow. I suppose that's the original purpose of having
# both parquet and csv. If so then I should also write GeoJSON when arrow is not
# installed. But you have to use arrow for newer versions of Vizgen in order to
# read the segmentations. You need GDAL to read GeoJSON.

test_that("readXenium XOA v1 flip image", {
    fp <- tempdir()
    fn <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
    sfe <- readXenium(fn, add_molecules = TRUE, flip = "image")
    # That things are aligned
    bbox_rg <- st_bbox(txSpots(sfe)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))

    img <- toExtImage(getImg(sfe), resolution = 1L)
    mask <- img > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe)))
    expect_true(mean(v$lyr.1) > 0.9)
    unlink(fn, recursive = TRUE)
})

# Read Xenium XOA v2===================
test_that("readXenium XOA v2, normal stuff", {
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    # Weirdly the first time I get the null pointer error
    sfe <- readXenium(fn, add_molecules = TRUE)
    # Basic stuff
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), c("centroids", "cellSeg", "nucSeg"))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(as.character(st_geometry_type(SpatialFeatureExperiment::centroids(sfe), FALSE)), "POINT")
    expect_equal(as.character(st_geometry_type(nucSeg(sfe), FALSE)), "MULTIPOLYGON")
    expect_equal(as.character(st_geometry_type(txSpots(sfe), FALSE)), "MULTIPOINT")
    expect_equal(imageIDs(sfe), "morphology_focus")
    expect_s4_class(getImg(sfe, image_id = "morphology_focus"), "BioFormatsImage")
    expect_equal(dim(getImg(sfe, image_id = "morphology_focus"))[["C"]], 4L)

    # That things are aligned
    bbox_rg <- st_bbox(txSpots(sfe)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))

    img <- toExtImage(getImg(sfe), resolution = 1L)
    mask <- img[,,1] > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe)))
    # About 1% of cells detected don't have nuclei here
    expect_true(mean(v$lyr.1, na.rm = TRUE) > 0.89)
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v2, somes images not found", {
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    unlink(file.path(fn, "morphology_focus"), recursive = TRUE)
    expect_warning(sfe <- readXenium(fn), "morphology_focus images not found")
    expect_true(isEmpty(imgData(sfe)))
    unlink(fn, recursive = TRUE)
})

test_that("readXenium XOA v2, use csv files", {
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    file.remove(list.files(fn, pattern = "*.parquet$", full.names = TRUE))
    expect_message(sfe <- readXenium(fn, add_molecules = TRUE),
                   ">>> Cell segmentations are found in `.csv` file")
    # Basic stuff
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), c("centroids", "cellSeg", "nucSeg"))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(as.character(st_geometry_type(SpatialFeatureExperiment::centroids(sfe), FALSE)), "POINT")
    expect_equal(as.character(st_geometry_type(nucSeg(sfe), FALSE)), "MULTIPOLYGON")
    expect_equal(as.character(st_geometry_type(txSpots(sfe), FALSE)), "MULTIPOINT")
    expect_equal(imageIDs(sfe), "morphology_focus")
    expect_s4_class(getImg(sfe, image_id = "morphology_focus"), "BioFormatsImage")
    expect_equal(dim(getImg(sfe, image_id = "morphology_focus"))[["C"]], 4L)

    # That things are aligned
    bbox_rg <- st_bbox(txSpots(sfe)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))

    img <- toExtImage(getImg(sfe), resolution = 1L)
    mask <- img[,,1] > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe)))
    # About 1% of cells detected don't have nuclei here
    expect_true(mean(v$lyr.1, na.rm = TRUE) > 0.89)
    unlink(fn, recursive = TRUE)
})

test_that("readXenium, flip image", {
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    sfe <- readXenium(fn, add_molecules = TRUE, flip = "image")
    sfe0 <- readXenium(fn)

    # That things are aligned
    bbox_rg <- st_bbox(txSpots(sfe)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))

    img <- toExtImage(getImg(sfe), resolution = 1L)
    mask <- img[,,1] > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe)))
    # About 1% of cells detected don't have nuclei here
    expect_true(mean(v$lyr.1, na.rm = TRUE) > 0.89)

    # That the image was actually flipped
    bfi <- getImg(sfe)
    expect_equal(transformation(bfi), list(name = "mirror", direction = "vertical"))
    unlink(fn, recursive = TRUE)
})

# Final cleanup in case failed test messed with cleanup
fp <- tempdir()
unlink(file.path(fp, "cosmx_test"), recursive = TRUE)
unlink(file.path(fp, "vizgen_test"), recursive = TRUE)
unlink(file.path(fp, "xenium_test"), recursive = TRUE)
