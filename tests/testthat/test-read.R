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

dir_use <- system.file("extdata/vizgen_cellbound", package = "SpatialFeatureExperiment")

test_that("readVizgen flip geometry, use cellpose", {
    file.copy(dir_use, ".", recursive = TRUE)
    sfe <- readVizgen("vizgen_cellbound", z = 3L, use_cellpose = TRUE,
                      flip = "geometry", min_area = 15)
    expect_equal(unit(sfe), "micron")
    expect_equal(imgData(sfe)$image_id,
                 paste0(c(paste0("Cellbound", 1:3), "DAPI", "PolyT"),
                       "_z3"))
    img <- imgRaster(getImg(sfe, image_id = "PolyT_z3"))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    v <- aggregate(v$mosaic_PolyT_z3, by = list(v$ID), FUN = sum)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$x, use = "complete.obs") > 0.4)
    # Make sure both segmentations and centroids are flipped
    hulls <- st_convex_hull(cellSeg(sfe))
    expect_true(all(vapply(seq_len(nrow(cg)), function(i) {
        st_covered_by(cg[i,], hulls[i,], sparse = FALSE)[1,1]
    }, FUN.VALUE = logical(1))))
    # Make sure that cells that are too small are removed
    cg <- cellSeg(sfe)
    areas <- st_area(cg)
    expect_true(all(areas > 15))
    unlink("vizgen", recursive = TRUE)
})

test_that("readVizgen write parquet file after reading hdf5", {
    file.copy(dir_use, ".", recursive = TRUE)
    file.remove(file.path("vizgen", "cellpose_micron_space.parquet"))
    sfe <- readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "PolyT")
    expect_true(file.exists(file.path("vizgen","hdf5s_micron_space.parquet")))
    # Second time reading
    expect_message(readVizgen("vizgen", z = 0L, image = "PolyT"),
                   "processed hdf5 files will be used")
    unlink("vizgen", recursive = TRUE)
})

test_that("readVizgen flip geometry, don't use cellpose", {
    file.copy(dir_use, ".", recursive = TRUE)
    sfe <- readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "PolyT",
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
    unlink("vizgen", recursive = TRUE)
})

test_that("readVizgen flip image", {
    file.copy(dir_use, ".", recursive = TRUE)
    sfe <- readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image")
    expect_equal(unit(sfe), "micron")
    img <- imgRaster(getImg(sfe))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mosaic_PolyT_z0) > 0.4)
    unlink("vizgen", recursive = TRUE)
})

test_that("readVizgen don't flip image when image is too large", {
    file.copy(dir_use, ".", recursive = TRUE)
    expect_error(readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "PolyT",
                            flip = "image", max_flip = "0.02 TB"),
                 "max_flip must be in either MB or GB")
    sfe <- readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image", max_flip = "0.02 MB")
    suppressWarnings(img_orig <- rast(file.path("vizgen", "images", "mosaic_PolyT_z0.tif")))
    img <- imgRaster(getImg(sfe))
    # Make sure image was not flipped
    expect_equal(terra::values(img), terra::values(img_orig))
    cg <- SpatialFeatureExperiment::centroids(sfe)
    v <- terra::extract(img, cg)
    nCounts <- Matrix::colSums(counts(sfe))
    expect_true(cor(nCounts, v$mosaic_PolyT_z0) > 0.4)
    unlink("vizgen", recursive = TRUE)
})

test_that("Don't flip image if it's GeoTIFF", {
    file.copy(dir_use, ".", recursive = TRUE)
    sfe <- readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "PolyT",
                      flip = "image")
    terra::writeRaster(imgRaster(getImg(sfe)),
                       filename = file.path("vizgen", "images", "mosaic_DAPI_z0.tif"),
                       overwrite = TRUE)
    sfe2 <- readVizgen("vizgen", z = 0L, use_cellpose = FALSE, image = "DAPI",
                       flip = "image")
    expect_equal(terra::values(imgRaster(getImg(sfe))), terra::values(imgRaster(getImg(sfe2))))
    file.remove(file.path("vizgen", "images", "mosaic_DAPI_z0.tif"))
    unlink("vizgen", recursive = TRUE)
})

test_that("Errors and warnings in readVizgen", {
    file.copy(dir_use, ".", recursive = TRUE)
    expect_warning(sfe <- readVizgen("vizgen", z = 0L, image = "DAPI", use_cellpose = FALSE),
                   "don't exist")
    expect_equal(nrow(imgData(sfe)), 0L)
    expect_error(readVizgen("vizgen", z = 7L, image = "PolyT", use_cellpose = FALSE),
                 "z must be beween 0 and 6")
    unlink("vizgen", recursive = TRUE)
})

# Make toy examples of multiple pieces
dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
parq <- sfarrow::st_read_parquet(file.path(dir_use, "cellpose_micron_space.parquet"))
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
    parq2 <- parq[1:4,]
    new_geo <- st_sfc(large_g, large_small, small_small, large_large)
    parq2$Geometry <- new_geo
    dir.create("multi")
    dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
    file.copy(dir_use, "multi", recursive = TRUE)

    dir_use <- file.path("multi", "vizgen")
    file.remove(file.path(dir_use, "cellpose_micron_space.parquet"))
    suppressWarnings(sfarrow::st_write_parquet(parq2, file.path(dir_use, "cellpose_micron_space.parquet")))

    w <- capture_warnings(sfe <- readVizgen(dir_use, z = 0L, image = "PolyT"))
    expect_match(w, "Sanity check", all = FALSE)
    expect_match(w, "The largest piece is kept", all = FALSE)
    cg <- cellSeg(sfe)
    expect_equal(st_geometry_type(cg, by_geometry = "FALSE") |> as.character(), "POLYGON")
    expect_equal(colnames(sfe), parq2$EntityID[c(1,2,4)])
    areas <- st_area(cg)
    expect_true(all(vapply(areas, all.equal, target = st_area(large_g),
                           FUN.VALUE = logical(1))))
    unlink("multi", recursive = TRUE)
})

test_that("No polygons left", {
    # Like they're all too small, or when the polygon file is empty, unlikely
    # but otherwise we get a mysterious error
    parq2 <- parq[1:2,]
    small <- st_polygon(small)
    new_geo <- st_sfc(small, small)
    parq2$Geometry <- new_geo
    dir.create("small")
    dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
    file.copy(dir_use, "small", recursive = TRUE)

    dir_use <- file.path("small", "vizgen")
    file.remove(file.path(dir_use, "cellpose_micron_space.parquet"))
    suppressWarnings(sfarrow::st_write_parquet(parq2, file.path(dir_use, "cellpose_micron_space.parquet")))

    expect_error(readVizgen(dir_use, z = 0L, image = "PolyT"),
                 "No polygons left after filtering")
    unlink("small", recursive = TRUE)
})

test_that("Image can be named polyT in older version", {
    dir.create("polyT")
    dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
    file.copy(dir_use, "polyT", recursive = TRUE)

    dir_use <- file.path("polyT", "vizgen")
    file.rename(file.path(dir_use, "images", "mosaic_PolyT_z0.tif"),
                file.path(dir_use, "images", "mosaic_polyT_z0.tif"))
    expect_no_warning(readVizgen(dir_use, z = 0L, image = "PolyT"))
    unlink("polyT", recursive = TRUE)
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
# Error messages when coordinate columns are not found
test_that("Read MERFISH transcript spots into rowGeometries", {
    dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
    dir.create("test_spots")
    file.copy(list.files(dir_use, full.names = TRUE), "test_spots", recursive = TRUE)
    sfe <- readVizgen("test_spots", z = 0L, image = "PolyT", add_molecules = TRUE)
    expect_equal(rowGeometryNames(sfe), "txSpots")
    rg <- txSpots(sfe)
    expect_equal(as.character(st_geometry_type(rg, FALSE)), "MULTIPOINT")
    unlink("test_spots", recursive = TRUE)
})

test_that("Format MERFISH transcript spots for colGeometries", {
    # Error when column in cell_col is absent
    dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
    expect_error(formatTxSpots(file.path(dir_use, "detected_transcripts.csv"),
                               dest = "colGeometry"),
                 "file_out must be specified")
    expect_error(formatTxSpots(file.path(dir_use, "detected_transcripts.csv"),
                               dest = "colGeometry", file_out = "vizgen/tx_in_cells"),
                 "Column indicating cell ID not found.")

    dir.create("cg_vizgen")
    file.copy(list.files(dir_use, full.names = TRUE), "cg_vizgen", recursive = TRUE)
    df <- read.csv(file.path("cg_vizgen", "detected_transcripts.csv"),
                   header = TRUE)
    sfe <- readVizgen(dir_use, z = 0L, image = "PolyT")
    df$cell_id <- sample(colnames(sfe), nrow(df), replace = TRUE)
    rownames(df) <- df$X
    df$X <- NULL
    write.csv(df, file = file.path("cg_vizgen", "detected_transcripts.csv"),
              row.names = TRUE, quote = FALSE)

    # First run
    cg <- formatTxSpots(file.path("cg_vizgen", "detected_transcripts.csv"),
                        dest = "colGeometry",
                        file_out = file.path("cg_vizgen", "tx_in_cells"),
                        z = 0L)
    dir_check <- file.path("cg_vizgen", "tx_in_cells")
    expect_equal(cg, dir_check)
    expect_true(dir.exists(dir_check))
    fns_expect <- paste0(unique(df$gene[df$global_z == 0L]), "_spots.parquet")
    fns <- list.files(dir_check)
    expect_setequal(fns, fns_expect)

    # Check contents
    fn <- file.path(dir_check, fns_expect[1])
    g <- sfarrow::st_read_parquet(fn)
    expect_equal(as.character(st_geometry_type(g, FALSE)), "MULTIPOINT")

    # Second run
    time_note <- Sys.time() # Check the files weren't written again
    cg <- formatTxSpots(file.path("cg_vizgen", "detected_transcripts.csv"),
                        dest = "colGeometry",
                        file_out = file.path("cg_vizgen", "tx_in_cells"))
    expect_equal(cg, normalizePath(dir_check))
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    unlink("cg_vizgen", recursive = TRUE)
})

# Error message when multiple parquet files are found in readVizgen

test_that("Error messages in formatTxSpots", {
    dir_use <- system.file("extdata/vizgen", package = "SpatialFeatureExperiment")
    file <- file.path(dir_use, "detected_transcripts.csv")
    expect_error(formatTxSpots(file, dest = "colGeometry", file_out = NULL),
                 "file_out must be specified")
    file2 <- system.file("extdata/sfe_visium.rds", package = "SpatialFeatureExperiment")
    expect_error(formatTxSpots(file2),
                 "The file must be one of csv, tsv, txt, or parquet")
    expect_error(formatTxSpots(file, z = "foo"),
                 "z must either be numeric")
    expect_error(formatTxSpots(file, spatialCoordsNames = c("foo", "bar")),
                 "foo, bar not found")
    expect_error(formatTxSpots(file, z = 8L),
                 "z plane\\(s\\) specified not found")
    expect_error(formatTxSpots(file, dest = "colGeometry", cell_col = "foo",
                               file_out = "bar"),
                 "Column indicating cell ID not found")
})

test_that("readCosMX, not reading transcript spots", {
    dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
    sfe <- readCosMX(dir_use, z = 1L)
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), c("centroids", "cellSeg"))
    expect_equal(st_geometry_type(cellSeg(sfe), by_geometry = FALSE) |> as.character(),
                 "POLYGON")
    expect_equal(dim(sfe), c(960, 27))
})

test_that("readCosMX, reading spots, 1 z-plane", {
    dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
    dir.create("cosmx")
    file.copy(list.files(dir_use, full.names = TRUE), "cosmx")

    sfe <- readCosMX("cosmx", z = 1L, add_molecules = TRUE)
    fn <- file.path("cosmx", "tx_spots.parquet")
    expect_true(file.exists(fn))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(st_geometry_type(txSpots(sfe), FALSE) |> as.character(),
                 "MULTIPOINT")
    expect_null(st_z_range(txSpots(sfe)))

    # Reloading the second time reading
    time_note <- Sys.time() # Check the files weren't written again
    sfe <- readCosMX("cosmx", z = 1L, add_molecules = TRUE)
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(st_geometry_type(txSpots(sfe), FALSE) |> as.character(),
                 "MULTIPOINT")
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    unlink("cosmx", recursive = TRUE)
})

test_that("readCosMX, 2 z-planes, split z, not splitting by cell compartments", {
    dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
    dir.create("cosmx")
    file.copy(list.files(dir_use, full.names = TRUE), "cosmx")

    sfe <- readCosMX("cosmx", z = "all", add_molecules = TRUE) # default split z
    d <- file.path("cosmx", "tx_spots")
    expect_true(dir.exists(d))
    fns_expect <- paste0("tx_spots_z", 0:1, ".parquet")
    expect_equal(list.files(d), fns_expect)
    expect_equal(rowGeometryNames(sfe), paste0("tx_spots_z", 0:1))
    expect_null(st_z_range(rowGeometry(sfe, "tx_spots_z0")))

    # Reloading the second time reading
    # Both z-planes
    time_note <- Sys.time()
    sfe <- readCosMX("cosmx", z = "all", add_molecules = TRUE)
    expect_null(st_z_range(rowGeometry(sfe, "tx_spots_z0")))
    fn <- file.path(d, "tx_spots_z0.parquet")
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    # Only read one of the z-planes
    sfe <- readCosMX("cosmx", z = 0L, add_molecules = TRUE)
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    unlink("cosmx", recursive = TRUE)
})

test_that("readCosMX, don't split z, don't split cell compartments", {
    dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
    dir.create("cosmx")
    file.copy(list.files(dir_use, full.names = TRUE), "cosmx")

    sfe <- readCosMX("cosmx", z = "all", add_molecules = TRUE,
                     z_option = "3d")
    fn <- file.path("cosmx", "tx_spots.parquet")
    expect_true(file.exists(fn))
    expect_equal(rowGeometryNames(sfe), "txSpots")
    expect_equal(unclass(st_z_range(txSpots(sfe))), c(zmin = 0, zmax = 1),
                 ignore_attr = "crs")

    unlink("cosmx", recursive = TRUE)
})


test_that("readCosMX, split z, split cell compartments", {
    dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
    dir.create("cosmx")
    file.copy(list.files(dir_use, full.names = TRUE), "cosmx")

    sfe <- readCosMX("cosmx", z = "all", add_molecules = TRUE,
                     split_cell_comps = TRUE)
    d <- file.path("cosmx", "tx_spots")
    expect_true(dir.exists(d))
    comps <- c("Nuclear", "None", "Membrane", "Cytoplasm")
    combs <- expand.grid(compartment = comps, z = 0:1, stringsAsFactors = FALSE)
    rgns <- paste0(combs$compartment, "_z", combs$z)
    fns <- paste0(rgns, ".parquet")

    expect_setequal(rowGeometryNames(sfe), rgns)
    expect_setequal(list.files(d), fns)

    # Reloading the second time reading
    time_note <- Sys.time()
    sfe <- readCosMX("cosmx", z = "all", add_molecules = TRUE,
                     split_cell_comps = TRUE)
    time_check <- file.info(file.path(d, fns[1]))$mtime
    expect_true(time_note > time_check)
    expect_setequal(rowGeometryNames(sfe), rgns)

    # Only read one of the z-planes
    sfe <- readCosMX("cosmx", z = 0L, add_molecules = TRUE,
                     split_cell_comps = TRUE)
    time_check <- file.info(file.path(d, fns[1]))$mtime
    expect_true(time_note > time_check)
    expect_setequal(rowGeometryNames(sfe), rgns[1:4])

    unlink("cosmx", recursive = TRUE)
})

test_that("Format CosMX spots for colGeometry, multiple z-planes", {
    dir_use <- system.file("extdata/cosmx", package = "SpatialFeatureExperiment")
    dir.create("cosmx")
    file.copy(list.files(dir_use, full.names = TRUE), "cosmx")

    cg <- formatTxSpots(file.path("cosmx", "Run5642_S3_Quarter_tx_file.csv"),
                        dest = "colGeometry", z = "all",
                        cell_col = c("cell_ID", "fov"),
                        gene_col = "target", not_in_cell_id = "0",
                        spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
                        file_out = file.path("cosmx", "tx_spots"))
    expect_equal(cg, file.path("cosmx", "tx_spots"))
    # Oh, great, there's Bex1/2, illegal file name. No wonder people don't like CosMX
    df <- data.table::fread(file.path("cosmx", "Run5642_S3_Quarter_exprMat_file.csv"))
    genes <- names(df)[-c(1:2)]
    combs <- expand.grid(gene = genes, z = 0:1, stringsAsFactors = FALSE)
    fns_expect <- paste0(gsub("/", ".", combs$gene), "_spots_z", combs$z, ".parquet")
    fns <- list.files(file.path("cosmx", "tx_spots"))
    # Not all genes have spots in this downsampled toy dataset
    expect_true(all(fns %in% fns_expect))

    fn <- file.path("cosmx", "tx_spots", fns[1])
    g <- sfarrow::st_read_parquet(fn)
    expect_equal(st_geometry_type(g, FALSE) |> as.character(), "MULTIPOINT")
    unlink("cosmx", recursive = TRUE)
})

unlink("cosmx", recursive = TRUE)
unlink("vizgen", recursive = TRUE)
