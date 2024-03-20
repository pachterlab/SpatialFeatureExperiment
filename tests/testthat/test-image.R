library(SFEData)
library(SpatialExperiment)
library(SingleCellExperiment)
library(terra)
library(RBioFormats)
library(EBImage)
library(sf)
library(xml2)

# SFE==================
sfe <- McKellarMuscleData("small")
img_path <- system.file(file.path("extdata", "sample01", "outs", "spatial",
                                  "tissue_lowres_image.png"),
                        package = "SpatialFeatureExperiment")
# TODO: Change later, copy entire directory in each test & delete afterwards
xenium_dir <- "~/SpatialFeatureExperiment/xenium_lr"
xenium_fn <- "~/SpatialFeatureExperiment/xenium_lr/morphology_mip.ome.tif"

test_that("addImg, SpatRasterImage", {
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        extent = "foo"),
                 "extent must be a numeric vector.")
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        extent = c(0, 14, 0, 23),
                        "must be of the set 'xmin', 'xmax', 'ymin', and 'ymax'"))
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        scale_fct = -1),
                 "scale_fct must be a positive number.")
    expect_error(addImg(sfe, "foo/bar", sample_id = "Vis5A", image_id = "ome"),
                 "img is not a valid file path.")
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    expect_s4_class(getImg(sfe), "SpatRasterImage")
    expect_error(addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                        scale_fct = 0.023),
                 "already contains an entry with")
})

test_that("addImg, BioFormatsImage", {
    skip_on_bioc()
    library(RBioFormats)
    # Weird, would always fail the first time
    try(addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome"))
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    expect_s4_class(getImg(sfe, image_id = "ome"), "BioFormatsImage")
})

test_that("addImg, EBImage", {
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
    skip_on_bioc()
    library(RBioFormats)
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    img <- getImg(sfe, image_id = "ome")
    sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "ome", resolution = 1L)
    img_t2 <- getImg(sfe, sample_id = "Vis5A", image_id = "ome")
    expect_equal(dim(imgRaster(img, resolution = 1L))[1:2], dim(imgRaster(img_t2))[2:1])
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
    skip_on_bioc()
    library(RBioFormats)
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    img <- getImg(sfe, image_id = "ome")
    img_m <- mirrorImg(img, resolution = 1L)
    mat1 <- as.array(imgRaster(img, resolution = 1L)) |> t()
    # SFE method
    sfe <- mirrorImg(sfe, sample_id = "Vis5A", image_id = "ome",
                     resolution = 1L)
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
    center <- bbox_center(ext_new)
    expect_equal(center, bbox_center(ext_old))
    x_dist <- (ext_old["xmax"] - ext_old["xmin"])/2
    y_dist <- (ext_old["ymax"] - ext_old["ymin"])/2
    expect_equal(ext_new[["xmin"]], unname(center[1] - y_dist))
    expect_equal(ext_new[["xmax"]], unname(center[1] + y_dist))
    expect_equal(ext_new[["ymin"]], unname(center[2] - x_dist))
    expect_equal(ext_new[["ymax"]], unname(center[2] + x_dist))
})

test_that("Rotate method, SFE, BioFormatsImage", {
    skip_on_bioc()
    library(RBioFormats)
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    bfi <- getImg(sfe, image_id = "ome")
    sfe <- rotateImg(sfe, sample_id = "Vis5A", image_id = "ome",
                     degrees = 90, resolution = 1L)
    ebi <- getImg(sfe, image_id = "ome")
    ext_new <- ext(ebi)
    ext_old <- ext(bfi)
    center <- bbox_center(ext_new)
    expect_equal(center, bbox_center(ext_old))
    x_dist <- (ext_old["xmax"] - ext_old["xmin"])/2
    y_dist <- (ext_old["ymax"] - ext_old["ymin"])/2
    expect_equal(ext_new[["xmin"]], unname(center[1] - y_dist))
    expect_equal(ext_new[["xmax"]], unname(center[1] + y_dist))
    expect_equal(ext_new[["ymin"]], unname(center[2] - x_dist))
    expect_equal(ext_new[["ymax"]], unname(center[2] + x_dist))
})

test_that("translateImg method, SFE", {
    v <- c(135, 246)
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    spi <- getImg(sfe, image_id = "lowres")
    sfe <- translateImg(sfe, sample_id = "Vis5A", image_id = "lowres",
                        v = v)
    spi_tr <- getImg(sfe, image_id = "lowres")
    ext_new <- ext(spi_tr)
    ext_old <- ext(spi)
    expect_equal(ext_new[c("xmin", "xmax")], ext_old[c("xmin", "xmax")] + v[1])
    expect_equal(ext_new[c("ymin", "ymax")], ext_old[c("ymin", "ymax")] + v[2])
})

# SpatRasterImage-------------
test_that("ext when the image is wrapped", {
    fn <- system.file(file.path("extdata", "sample01", "outs", "spatial",
                                "tissue_hires_image.png"),
                      package = "SpatialFeatureExperiment")
    ext_use <- c(xmin = 0, xmax = 45, ymin = 0, ymax = 77)
    suppressWarnings(r <- rast(fn))
    rr <- wrap(r)
    spi <- SpatRasterImage(rr)
    ext1 <- ext(spi)
    expect_equal(ext1, ext_use)
})

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
    center <- bbox_center(ext_new)
    expect_equal(center, bbox_center(ext_old))
    x_dist <- (ext_old["xmax"] - ext_old["xmin"])/2
    y_dist <- (ext_old["ymax"] - ext_old["ymin"])/2
    expect_equal(ext_new[["xmin"]], unname(center[1] - y_dist))
    expect_equal(ext_new[["xmax"]], unname(center[1] + y_dist))
    expect_equal(ext_new[["ymin"]], unname(center[2] - x_dist))
    expect_equal(ext_new[["ymax"]], unname(center[2] + x_dist))
})

test_that("translateImg, SpatRasterImage method", {
    img <- SpatRasterImage(suppressWarnings(rast(img_path)))
    ext_old <- ext(img)
    v <- c(135, 246)
    img_tr <- translateImg(img, v)
    ext_new <- ext(img_tr)
    expect_equal(ext_new[c("xmin", "xmax")], ext_old[c("xmin", "xmax")] + v[1])
    expect_equal(ext_new[c("ymin", "ymax")], ext_old[c("ymin", "ymax")] + v[2])
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
    expect_equal(imgRaster(spi) |> terra::as.array(),
                 imgRaster(ebi) |> as.array() |> aperm(c(2,1,3)))
})
fn <- system.file(file.path("extdata", "vizgen_cellbound", "images",
                            "mosaic_Cellbound1_z3.tif"),
                  package = "SpatialFeatureExperiment")
test_that("Convert SpatRasterImage to EBImage, grayscale", {
    suppressWarnings(spi <- SpatRasterImage(rast(fn)))
    ebi <- toEBImage(spi)
    expect_equal(ext(spi), ext(ebi))
    m1 <- imgRaster(spi) |> terra::as.array()
    m1 <- m1[,,1]
    expect_equal(m1, imgRaster(ebi) |> as.array() |> t())
})

# BioFormatsImage=================
xml_meta <- read.omexml(xenium_fn) |> read_xml() |> as_list()
psx <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeX") |> as.numeric()
psy <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeY") |> as.numeric()
meta <- read.metadata(xenium_fn) |> coreMetadata(1)
sizeX_full <- meta$sizeX
sizeY_full <- meta$sizeY
test_that("BioFormatsImage constructor", {
    skip_on_bioc()
    # Extent inferred from metadata
    bfi <- BioFormatsImage(xenium_fn)
    expect_s4_class(bfi, "BioFormatsImage")
    expect_true(isFull(bfi))
    expect_equal(imgSource(bfi), xenium_fn)
    ext_bfi <- ext(bfi)
    expect_setequal(names(ext_bfi), c("xmin", "xmax", "ymin", "ymax"))
    ext_expect <- c(xmin = 0, ymin = 0,
                    xmax = sizeX_full*psx, ymax = sizeY_full*psy)
    expect_equal(ext_bfi, ext_expect[names(ext_bfi)])
})

test_that("dim for BioFormatsImage", {
    bfi <- BioFormatsImage(xenium_fn)
    d <- dim(bfi)
    expect_equal(d, c(sizeX_full, sizeY_full))
})

test_that("Errors when constructing BioFormatsImage", {
    skip_on_bioc()
    library(RBioFormats)
    # Invalid bbox
    bbox_use <- c(xmin = 100, xmax = 0, ymin = 0, ymax = 102)
    expect_error(BioFormatsImage(xenium_fn, ext = bbox_use, isFull = FALSE),
                 "Min limit is greater than max limit")
    # isFull
    expect_error(BioFormatsImage(xenium_fn, ext = bbox_use, isFull = NA),
                 "isFull must be either TRUE or FALSE, not NA")
})

test_that("Convert BioFormatsImage to EBImage, full extent", {
    skip_on_bioc()
    library(RBioFormats)
    ext_expect <- c(xmin = 0, ymin = 0, xmax = sizeX_full*psx, ymax = sizeY_full*psy)
    bfi <- BioFormatsImage(xenium_fn)
    expect_error(toEBImage(bfi, resolution = 10L),
                 "Resolution subscript out of bound")
    ebi <- toEBImage(bfi, resolution = 1L)
    expect_s4_class(ebi, "EBImage")
    expect_equal(ext(ebi), ext_expect[c("xmin", "xmax", "ymin", "ymax")])
    dim_img <- dim(imgRaster(ebi))
    expect_equal(dim_img, c(sizeX_full, sizeY_full))
})
ext_use <- c(xmin = 1000, xmax = 2000, ymin = 600, ymax = 1600)
test_that("Convert BioFormatsImage to EBImage, not full extent", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    ebi <- toEBImage(bfi, resolution = 1L)
    expect_s4_class(ebi, "EBImage")
    # Extent from lower res must be slightly larger as I include the boundary pixels
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_ebi <- st_bbox(ext(ebi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_ebi, sparse = FALSE))
    expect_true(st_area(ext_sf_ebi) / st_area(ext_sf_bfi) < 1.005)
    # pixel range
    dim_expect <- c(588, 589)
    expect_equal(dim(imgRaster(ebi)), dim_expect)
    # Check content, rather crude, check that it includes that big bright patch
    expect_true(sum(imgRaster(ebi) > 1e4) > 400)
    expect_true(all(imgRaster(ebi)[199:222,440:448]))
})

test_that("Convert BioFormatsImage to EBImage, not resolution 1", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    ebi <- toEBImage(bfi, resolution = 2L)
    expect_s4_class(ebi, "EBImage")
    # Extent from lower res must be slightly larger as I include the boundary pixels
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_ebi <- st_bbox(ext(ebi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_ebi, sparse = FALSE))
    expect_true(st_area(ext_sf_ebi) / st_area(ext_sf_bfi) < 1.01)
    # pixel range
    dim_expect <- c(196, 196)
    expect_equal(dim(imgRaster(ebi)), dim_expect)
    # Check content, rather crude, check that it includes that big bright patch
    expect_true(sum(imgRaster(ebi) > 1e4) > 50)
    expect_true(all(imgRaster(ebi)[67:75,148:151]))
})

test_that("When physical pixel size is absent from metadata", {
    skip_on_bioc()
    fn <- system.file(file.path("extdata", "xenium", "morphology_focus.ome.tif"),
                      package = "SpatialFeatureExperiment")
    ext_use <- c(xmin = 0, xmax = 237, ymin = 0, ymax = 237)
    expect_warning(bfi <- BioFormatsImage(fn),
                   "Physical pixel size absent from image metadata")
    expect_equal(ext(bfi), ext_use)
})

test_that("Ignore resolution in toEBImage when there's only 1 resolution", {
    skip_on_bioc()
    # TODO: change file path after I make the toy dataset
    fn <- system.file(file.path("extdata", "xenium", "morphology_focus.ome.tif"),
                      package = "SpatialFeatureExperiment")
    ext_use <- c(xmin = 0, xmax = 237, ymin = 0, ymax = 237)
    suppressWarnings(bfi <- BioFormatsImage(fn))
    expect_warning(ebi <- toEBImage(bfi, resolution = 4L),
                   "Physical pixel size absent from image metadata")
    expect_equal(ext(ebi), ext_use)
})

test_that("Convert BioFormatsImage to SpatRasterImage", {
    skip_on_bioc()
    # TODO: After making smaller xenium subset, copy the whole directory and delete after test
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    expect_message(spi <- toSpatRasterImage(bfi, resolution = 1L), "non OME")
    fn <- file.path(xenium_dir, "morphology_mip_res1.tiff")
    expect_true(file.exists(fn))
    expect_s4_class(spi, "SpatRasterImage")
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_spi <- st_bbox(ext(spi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_spi, sparse = FALSE))
    expect_true(st_area(ext_sf_spi) / st_area(ext_sf_bfi) < 1.005)
    unlink(fn)
})
test_that("Convert BioFormatsImage to SpatRasterImage not saving geotiff", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    expect_no_message(spi <- toSpatRasterImage(bfi, save_geotiff = FALSE,
                                               resolution = 1L))
    fn <- file.path(xenium_dir, "morphology_mip_res1.tiff")
    expect_false(file.exists(fn))
    expect_s4_class(spi, "SpatRasterImage")
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_spi <- st_bbox(ext(spi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_spi, sparse = FALSE))
    expect_true(st_area(ext_sf_spi) / st_area(ext_sf_bfi) < 1.005)
})

ext_use2 <- c(xmin = 200, xmax = 1700, ymin = 200, ymax = 1200)
test_that("transpose, check ext", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- transposeImg(bfi, resolution = 1L)
    expect_s4_class(ebi, "EBImage")
    ext_ebi <- ext(ebi)
    expect_true(ext_ebi["ymax"] > ext_ebi["xmax"])
    expect_equal(imgRaster(ebi), toEBImage(bfi, resolution = 1L) |> imgRaster() |> EBImage::transpose())
})

test_that("mirror (EBI behind the scene), vertical", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- mirrorImg(bfi, direction = "vertical", resolution = 1L)
    expect_s4_class(ebi, "EBImage")
    ext_ebi <- ext(ebi)
    ebi_orig <- toEBImage(bfi, resolution = 1L)
    expect_equal(ext_ebi, ext(ebi_orig))
    expect_equal(imgRaster(ebi), EBImage::flip(imgRaster(ebi_orig)))
})

test_that("mirror (EBI behind the scene), horizontal", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- mirrorImg(bfi, direction = "horizontal", resolution = 1L)
    expect_s4_class(ebi, "EBImage")
    ext_ebi <- ext(ebi)
    ebi_orig <- toEBImage(bfi, resolution = 1L)
    expect_equal(ext_ebi, ext(ebi_orig))
    expect_equal(imgRaster(ebi), EBImage::flop(imgRaster(ebi_orig)))
})

test_that("rotate (EBI behind the scene)", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    expect_error(rotateImg(bfi, 45), "degrees%%90 == 0 is not TRUE")
    ebi <- rotateImg(bfi, 90, resolution = 1L)
    ext_ebi <- ext(ebi)
    ebi_orig <- toEBImage(bfi, resolution = 1L)
    expect_equal(imgRaster(ebi), EBImage::rotate(imgRaster(ebi_orig), 90))
    ext_new <- ext(ebi)
    ext_old <- ext(ebi_orig)
    center <- bbox_center(ext_new)
    expect_equal(center, bbox_center(ext_old))
    x_dist <- (ext_old["xmax"] - ext_old["xmin"])/2
    y_dist <- (ext_old["ymax"] - ext_old["ymin"])/2
    bbox_exp <- c(center[1] - y_dist, center[1] + y_dist,
                  center[2] - x_dist, center[2] + x_dist)
    names(bbox_exp) <- c("xmin", "xmax", "ymin", "ymax")
    new_sf <- st_bbox(ext_new) |> st_as_sfc()
    exp_sf <- st_bbox(bbox_exp) |> st_as_sfc()
    expect_true(st_covered_by(new_sf, exp_sf, sparse = FALSE))
    expect_equal(st_area(new_sf), st_area(exp_sf))
})

test_that("translateImg, BioFormatsImage method", {
    skip_on_bioc()
    bfi <- BioFormatsImage(xenium_fn)
    v <- c(135, 246)
    ext_old <- ext(bfi)
    bfi_tr <- translateImg(bfi, v)
    ext_new <- ext(bfi_tr)
    expect_equal(ext_new[c("xmin", "xmax")], ext_old[c("xmin", "xmax")] + v[1])
    expect_equal(ext_new[c("ymin", "ymax")], ext_old[c("ymin", "ymax")] + v[2])
    expect_equal(SpatialFeatureExperiment::origin(bfi_tr), v)
})

test_that("translateImg, EBImage method", {
    skip_on_bioc()
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- toEBImage(bfi, resolution = 1L)
    ext_old <- ext(ebi)
    v <- c(135, 246)
    ebi_tr <- translateImg(ebi, v)
    ext_new <- ext(ebi_tr)
    expect_equal(ext_new[c("xmin", "xmax")], ext_old[c("xmin", "xmax")] + v[1])
    expect_equal(ext_new[c("ymin", "ymax")], ext_old[c("ymin", "ymax")] + v[2])
})

test_that("Convert BioFormatsImage to EBImage after translation", {
    skip_on_bioc()
    library(RBioFormats)
    # Make sure that the right part of the image is loaded
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- toEBImage(bfi, resolution = 1L)
    v <- c(135, 246)
    bfi_tr <- translateImg(bfi, v)
    ebi_sub <- toEBImage(bfi_tr, resolution = 1L)
    ext_sub <- ext(ebi_sub)
    ext_exp <- ext(ebi)
    ext_exp[c("xmin", "xmax")] <- ext_exp[c("xmin", "xmax")] + v[1]
    ext_exp[c("ymin", "ymax")] <- ext_exp[c("ymin", "ymax")] + v[2]
    expect_equal(ext_sub, ext_exp)
    expect_equal(imgRaster(ebi), imgRaster(ebi_sub))
})

test_that("Crop SpatRasterImage", {
    suppressWarnings(spi <- SpatRasterImage(rast(fn)))
    bbox_use <- c(xmin = 100, xmax = 200, ymin = 100, ymax = 200)
    spi_sub <- cropImg(spi, bbox_use)
    expect_s4_class(spi_sub, "SpatRasterImage")
    expect_equal(ext(spi_sub), bbox_use)
})

test_that("Crop BioFormatsImage", {
    skip_on_bioc()
    library(RBioFormats)
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
    skip_on_bioc()
    library(RBioFormats)
    ebi <- BioFormatsImage(xenium_fn) |> toEBImage(resolution = 1L)
    ebi_sub <- cropImg(ebi, ext_use)
    expect_s4_class(ebi_sub, "EBImage")
    ext_sf1 <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf2 <- st_bbox(ext(ebi_sub)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf1, ext_sf2, sparse = FALSE))
    expect_true(st_area(ext_sf2) / st_area(ext_sf1) < 1.005)
    # pixel range
    dim_expect <- c(588, 589)
    expect_equal(dim(imgRaster(ebi_sub)), round(dim_expect))
})

test_that("EBImage constructor", {
    img <- readImage(system.file('images', 'nuclei.tif', package='EBImage'))
    expect_error(EBImage(img, ext = NULL),
                 "Extent must be specified")
})
