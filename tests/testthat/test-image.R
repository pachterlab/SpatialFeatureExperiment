library(SFEData)
library(SpatialExperiment)
library(SingleCellExperiment)
library(terra)
library(RBioFormats)
library(EBImage)
library(sf)

# SFE==================
sfe <- McKellarMuscleData("small")
img_path <- system.file(file.path("extdata", "sample01", "outs", "spatial",
                                  "tissue_lowres_image.png"),
                        package = "SpatialFeatureExperiment")
# TODO: Change later, copy entire directory in each test & delete afterwards
xenium_fn <- "xenium_toy/morphology_mip.ome.tif"

test_that("addImg", {
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        extent = "foo"),
                 "extent must be a numeric vector.")
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        extent = c(0, 14, 0, 23),
                        "must be of the set 'xmin', 'xmax', 'ymin', and 'ymax'"))
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        scale_fct = -1),
                 "scale_fct must be a positive number.")
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    expect_s4_class(getImg(sfe), "SpatRasterImage")
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        scale_fct = 0.023),
                 "already contains an entry with")

    # Add BioFormatsImage
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    expect_s4_class(getImg(sfe, image_id = "ome"), "BioFormatsImage")

    # Add EBImage
    img <- readImage(system.file('images', 'nuclei.tif', package='EBImage'))
    ext_img <- c(xmin = 0, xmax = dim(img)[1], ymin = 0, ymax = dim(img)[2])
    sfe <- addImg(sfe, img, sample_id = "Vis5A", image_id = "ebi", extent = ext_img)
    expect_s4_class(getImg(sfe, image_id = "ebi"), "EBImage")
})

test_that("transposeImg, SFE method, SpatRasterImage", {
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    img <- getImg(sfe)
    sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    img_t2 <- getImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    expect_equal(dim(imgRaster(img))[1:2], dim(imgRaster(img_t2))[2:1])
})

test_that("transposeImg, SFE method, BioFormatsImage", {
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    img <- getImg(sfe, image_id = "ome")
    sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "ome")
    img_t2 <- getImg(sfe, sample_id = "Vis5A", image_id = "ome")
    expect_equal(dim(imgRaster(img))[1:2], dim(imgRaster(img_t2))[2:1])
})

test_that("transposeImg, SFE method, EBImage", {
    img <- readImage(system.file('images', 'nuclei.tif', package='EBImage'))
    ext_img <- c(xmin = 0, xmax = dim(img)[1], ymin = 0, ymax = dim(img)[2])
    sfe <- addImg(sfe, img, sample_id = "Vis5A", image_id = "ebi", extent = ext_img)
    img <- getImg(sfe, image_id = "ebi")
    sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "ebi")
    img_t2 <- getImg(sfe, sample_id = "Vis5A", image_id = "ebi")
    expect_equal(dim(imgRaster(img))[1:2], dim(imgRaster(img_t2))[2:1])
})

test_that("mirrorImg, SFE method, SpatRasterImage", {
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    img <- getImg(sfe)
    img_m <- mirrorImg(img)
    mat1 <- terra::as.array(terra::mean(imgRaster(img)))[,,1]
    # SFE method
    sfe <- mirrorImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    img_m2 <- getImg(sfe)
    mat3 <- terra::as.array(terra::mean(imgRaster(img_m2)))[,,1]
    mat3_rev <- apply(mat3, 2, rev)
    expect_equal(mat1, mat3_rev)
})

test_that("mirrorImg, SFE method, BioFormatsImage", {
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    img <- getImg(sfe, image_id = "ome")
    img_m <- mirrorImg(img)
    mat1 <- as.array(imgRaster(img)) |> t()
    # SFE method
    sfe <- mirrorImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    img_m2 <- getImg(sfe)
    mat3 <- as.array(imgRaster(img_m)) |> t()
    mat3_rev <- apply(mat3, 2, rev)
    expect_equal(mat1, mat3_rev)
})

test_that("Rotate method, SFE, SpatRasterImage", {
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    spi <- getImg(sfe, image_id = "lowres")
    sfe <- rotateImg(sfe, sample_id = "Vis5A", image_id = "lowres",
                     degrees = 90)
    ebi <- getImg(sfe, image_id = "lowres")
    ext_new <- ext(ebi)
    ext_old <- ext(spi)
    expect_equal(ext_new[["ymax"]], ext_old[["xmax"]])
    expect_equal(ext_new[["xmax"]], ext_old[["ymax"]])
})

test_that("Rotate method, SFE, BioFormatsImage", {
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    bfi <- getImg(sfe, image_id = "ome")
    sfe <- rotateImg(sfe, sample_id = "Vis5A", image_id = "ome",
                     degrees = 90)
    ebi <- getImg(sfe, image_id = "ome")
    ext_new <- ext(ebi)
    ext_old <- ext(bfi)
    expect_equal(ext_new[["ymax"]], ext_old[["xmax"]])
    expect_equal(ext_new[["xmax"]], ext_old[["ymax"]])
})

# SpatRasterImage-------------
test_that("transposeImg, SpatRasterImage method", {
    img <- SpatRasterImage(suppressWarnings(rast(img_path)))
    # SpatRasterImage method
    img_t <- transposeImg(img)
    expect_s4_class(img_t, "SpatRasterImage")
    expect_equal(dim(imgRaster(img))[1:2], dim(imgRaster(img_t))[2:1])
})

test_that("mirrorImg, SpatRasterImage method", {
    img <- SpatRasterImage(suppressWarnings(rast(img_path)))
    img_m <- mirrorImg(img)
    mat1 <- terra::as.array(terra::mean(imgRaster(img)))[,,1]
    mat2 <- terra::as.array(terra::mean(imgRaster(img_m)))[,,1]
    mat2_rev <- apply(mat2, 2, rev)
    expect_equal(mat1, mat2_rev)
})

test_that("Rotate method for SpatRasterImage which converts to EBImage", {
    suppressWarnings(spi <- SpatRasterImage(rast(img_path)))
    expect_error(rotateImg(spi, 45), "degrees%%90 == 0 is not TRUE")
    ebi <- rotateImg(spi, 90)
    ext_new <- ext(ebi)
    ext_old <- ext(spi)
    expect_equal(ext_new[["ymax"]], ext_old[["xmax"]])
    expect_equal(ext_new[["xmax"]], ext_old[["ymax"]])
})

sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
              scale_fct = 0.023)
test_that("imgRaster, trivial", {
    img <- imgRaster(getImg(sfe))
    expect_s4_class(img, "SpatRaster")
})

test_that("imgRaster, loaded from RDS", {
    saveRDS(sfe, "baz.rds")
    sfe_read <- readRDS("baz.rds")
    img <- imgRaster(getImg(sfe))
    expect_s4_class(img, "SpatRaster")
    unlink("baz.rds")
})

test_that("imgSource, trivial", {
    expect_type(imgSource(getImg(sfe)), "character")
})

test_that("Convert SpatRasterImage to EBImage, RGB", {
    suppressWarnings(spi <- SpatRasterImage(rast(img_path)))
    ebi <- toEBImage(spi)
    expect_s4_class(ebi, "EBImage")
    expect_equal(ext(spi), ext(ebi))
    expect_equal(imgRaster(spi) |> as.array(),
                 imgRaster(ebi) |> as.array() |> aperm(c(2,1,3)))
})

test_that("Convert SpatRasterImage to EBImage, grayscale", {
    fn <- "inst/extdata/vizgen_cellbound/images/mosaic_Cellbound1_z3.tif"
    suppressWarnings(spi <- SpatRasterImage(rast(fn)))
    ebi <- toEBImage(spi)
    expect_equal(ext(spi), ext(ebi))
    m1 <- imgRaster(spi) |> as.array()
    m1 <- m1[,,1]
    expect_equal(m1, imgRaster(ebi) |> as.array() |> t())
})

# BioFormatsImage=================
test_that("BioFormatsImage constructor", {
    # Extent inferred from metadata
    bfi <- BioFormatsImage(xenium_fn)
    expect_s4_class(bfi, "BioFormatsImage")
    expect_true(isFull(bfi))
    expect_equal(imgSource(bfi), xenium_fn)
    ext_bfi <- ext(bfi)
    expect_setequal(names(ext_bfi), c("xmin", "xmax", "ymin", "ymax"))
    ext_expect <- c(xmin = 0, ymin = 0, xmax = 11454*0.2125, ymax = 10307*0.2125)
    expect_equal(ext_bfi, ext_expect)
})

test_that("Errors when constructing BioFormatsImage", {
    # Invalid bbox
    bbox_use <- c(xmin = 100, xmax = 0, ymin = 0, ymax = 102)
    expect_error(BioFormatsImage(xenium_fn, ext = bbox_use, isFull = FALSE),
                 "Min limit is greater than max limit")
    # isFull
    expect_error(BioFormatsImage(xenium_fn, ext = bbox_use, isFull = NA),
                 "isFull must be either TRUE or FALSE, not NA")
    # Trying to read image that doesn't have physical pixel size
    expect_error(BioFormatsImage("inst/extdata/vizgen_cellbound/images/mosaic_Cellbound1_z3.tif"),
                 "Physical pixel size absent from image metadata.")
})

sizeX_full <- 11454
sizeY_full <- 10307
psx <- psy <- 0.2125
sizeX4 <- 1431
sizeY4 <- 1288
test_that("Convert BioFormatsImage to EBImage, full extent", {
    ext_expect <- c(xmin = 0, ymin = 0, xmax = sizeX_full*psx, ymax = sizeY_full*psy)
    bfi <- BioFormatsImage(xenium_fn)
    ebi <- toEBImage(bfi, resolution = 4L)
    expect_s4_class(ebi, "EBImage")
    expect_equal(ext(ebi), ext_expect)
    dim_img <- dim(imgRaster(ebi))
    expect_equal(dim_img, c(sizeX4, sizeY4))
})
ext_use <- c(xmin = 200, xmax = 1200, ymin = 200, ymax = 1200)
test_that("Convert BioFormatsImage to EBImage, not full extent", {
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    ebi <- toEBImage(bfi, resolution = 4L)
    expect_s4_class(ebi, "EBImage")
    # Extent from lower res must be slightly larger as I include the boundary pixels
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_ebi <- st_bbox(ext(ebi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_ebi, sparse = FALSE))
    expect_true(st_area(ext_sf_ebi) / st_area(ext_sf_bfi) < 1.005)
    # pixel range
    dim_expect <- c(1000/(sizeX_full*psx)*sizeX4, 1000/(sizeY_full*psy)*sizeY4)
    expect_equal(dim(imgRaster(ebi)), round(dim_expect))
})

test_that("Convert BioFormatsImage to SpatRasterImage", {
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    expect_message(spi <- toSpatRasterImage(bfi, resolution = 4L), "Saving image")
    fn <- "xenium_toy/morphology_mip_res4.tif"
    expect_true(file.exists(fn))
    expect_s4_class(spi, "SpatRasterImage")
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_spi <- st_bbox(ext(spi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_spi, sparse = FALSE))
    expect_true(st_area(ext_sf_spi) / st_area(ext_sf_bfi) < 1.005)
})
ext_use2 <- c(xmin = 200, xmax = 1700, ymin = 200, ymax = 1200)
test_that("transpose, check ext", {
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- transposeImg(bfi)
    expect_s4_class(ebi, "EBImage")
    ext_ebi <- ext(ebi)
    expect_true(ext_ebi["ymax"] > ext_ebi["xmax"])
    expect_equal(imgRaster(ebi), toEBImage(bfi) |> imgRaster() |> EBImage::transpose())
})

test_that("mirror (EBI behind the scene), vertical", {
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- mirrorImg(bfi, direction = "vertical")
    expect_s4_class(ebi, "EBImage")
    ext_ebi <- ext(ebi)
    ebi_orig <- toEBImage(bfi)
    expect_equal(ext_ebi, ext(ebi_orig))
    expect_equal(imgRaster(ebi), EBImage::flip(imgRaster(ebi_orig)))
})

test_that("mirror (EBI behind the scene), horizontal", {
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- mirrorImg(bfi, direction = "horizontal")
    expect_s4_class(ebi, "EBImage")
    ext_ebi <- ext(ebi)
    ebi_orig <- toEBImage(bfi)
    expect_equal(ext_ebi, ext(ebi_orig))
    expect_equal(imgRaster(ebi), EBImage::flop(imgRaster(ebi_orig)))
})

test_that("rotate (EBI behind the scene)", {
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    expect_error(rotateImg(bfi, 45), "degrees%%90 == 0 is not TRUE")
    ebi <- rotateImg(bfi, 90)
    ext_ebi <- ext(ebi)
    ebi_orig <- toEBImage(bfi)
    expect_true(ext_ebi["ymax"] > ext_ebi["xmax"])
    expect_equal(imgRaster(ebi), EBImage::rotate(imgRaster(ebi_orig), 90))
})

test_that("Crop SpatRasterImage", {
    fn <- "inst/extdata/vizgen_cellbound/images/mosaic_Cellbound1_z3.tif"
    suppressWarnings(spi <- SpatRasterImage(rast(fn)))
    bbox_use <- c(xmin = 100, xmax = 200, ymin = 100, ymax = 200)
    spi_sub <- cropImg(spi, bbox_use)
    expect_s4_class(spi_sub, "SpatRasterImage")
    expect_equal(ext(spi_sub), bbox_use)
})

test_that("Crop BioFormatsImage", {
    bfi <- BioFormatsImage(xenium_fn)
    bfi_sub <- cropImg(bfi, ext_use)
    expect_s4_class(bfi_sub, "BioFormatsImage")
    expect_equal(ext(bfi_sub), ext_use)
    # When bbox is not fully contained inside the original extent
    ext_test <- c(xmin = 2000, xmax = 3000, ymin = 2000, ymax = 3000)
    bfi_sub2 <- cropImg(bfi, ext_test)
    ext_old <- ext(bfi)
    expect_equal(ext(bfi_sub2),
                 c(xmin = 2000, xmax = ext_old[["xmax"]], ymin = 2000, ymax = ext_old[["ymax"]]))
})

test_that("Crop EBImage", {
    ebi <- BioFormatsImage(xenium_fn) |> toEBImage(resolution = 4L)
    ebi_sub <- cropImg(ebi, ext_use)
    expect_s4_class(ebi_sub, "EBImage")
    ext_sf1 <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf2 <- st_bbox(ext(ebi_sub)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf1, ext_sf2, sparse = FALSE))
    expect_true(st_area(ext_sf2) / st_area(ext_sf1) < 1.005)
    # pixel range
    dim_expect <- c(1000/(sizeX_full*psx)*sizeX4, 1000/(sizeY_full*psy)*sizeY4)
    expect_equal(dim(imgRaster(ebi_sub)), round(dim_expect))
})

test_that("EBImage constructor", {
    img <- readImage(system.file('images', 'nuclei.tif', package='EBImage'))
    expect_error(EBImage(img, ext = NULL),
                 "Extent must be specified")
})
