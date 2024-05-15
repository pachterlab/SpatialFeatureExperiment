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
fp <- tempdir()
xenium_dir <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
xenium_fn <- file.path(xenium_dir, "morphology_mip.ome.tif")

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
    library(RBioFormats)
    # Weird, would always fail the first time
    try(addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome"))
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    expect_s4_class(getImg(sfe, image_id = "ome"), "BioFormatsImage")
})

test_that("addImg, ExtImage", {
    img <- readImage(system.file('images', 'nuclei.tif', package="EBImage"))
    ext_img <- c(xmin = 0, xmax = dim(img)[1], ymin = 0, ymax = dim(img)[2])
    sfe <- addImg(sfe, img, sample_id = "Vis5A", image_id = "ebi", extent = ext_img)
    expect_s4_class(getImg(sfe, image_id = "ebi"), "ExtImage")
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
    library(RBioFormats)
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    img <- getImg(sfe, image_id = "ome")
    sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "ome")
    img_t2 <- getImg(sfe, sample_id = "Vis5A", image_id = "ome")
    expect_equal(unname(ext(img)), unname(ext(img_t2)[c("ymin", "ymax", "xmin", "xmax")]))
})

test_that("transposeImg, SFE method, ExtImage", {
    img <- readImage(system.file('images', 'nuclei.tif', package="EBImage"))
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
    library(RBioFormats)
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    img <- getImg(sfe, image_id = "ome")
    img_m <- mirrorImg(img)
    # SFE method
    sfe <- mirrorImg(sfe, sample_id = "Vis5A", image_id = "ome")
    img_m2 <- getImg(sfe, image_id = "ome")
    expect_equal(img_m, img_m2)
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
    library(RBioFormats)
    sfe <- addImg(sfe, xenium_fn, sample_id = "Vis5A", image_id = "ome")
    bfi <- getImg(sfe, image_id = "ome")
    sfe <- rotateImg(sfe, sample_id = "Vis5A", image_id = "ome",
                     degrees = 90)
    bfi2 <- getImg(sfe, image_id = "ome")
    ext_new <- ext(bfi2)
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
test_that("show for SpatRasterImage", {
    img <- SpatRasterImage(suppressWarnings(rast(img_path)))
    expect_output(show(img), "14 x 23 x 3 \\(width x height x channels\\) SpatRasterImage")
})

test_that("transposeImg, SpatRasterImage method", {
    img <- SpatRasterImage(suppressWarnings(rast(img_path)))
    # SpatRasterImage method
    img_t <- transposeImg(img)
    expect_s4_class(img_t, "SpatRasterImage")
    expect_equal(dim(imgRaster(img))[1:2], dim(imgRaster(img_t))[2:1])

    # Use maxcell
    img_t2 <- transposeImg(img, maxcell = 100)
    d <- dim(img_t2)
    expect_equal(d[1], 8)
    expect_equal(d[2], 13)

    # Save to file
    img_t3 <- transposeImg(img, filename = "foo.tif")
    expect_equal(dim(img)[1:2], dim(img_t3)[2:1])
    expect_equal(imgSource(img_t3), "foo.tif")
    expect_true(file.exists("foo.tif"))
    file.remove("foo.tif")
    file.remove("foo.tif.aux.xml")
})

test_that("mirrorImg, SpatRasterImage method", {
    img <- SpatRasterImage(suppressWarnings(rast(img_path)))
    img_m <- mirrorImg(img)
    mat1 <- terra::as.array(terra::mean(imgRaster(img)))[,,1]
    mat2 <- terra::as.array(terra::mean(imgRaster(img_m)))[,,1]
    mat2_rev <- apply(mat2, 2, rev)
    expect_equal(mat1, mat2_rev)

    # Use maxcell
    img_m2 <- mirrorImg(img, maxcell = 100)
    expect_equal(dim(img_m2), c(13,8,3))
    # check content
    expect_equal(imgRaster(img_m2)[10,3,2][[1]], 190)

    # Use filename
    img_m3 <- mirrorImg(img, filename = "foo.tif")
    expect_equal(imgSource(img_m3), "foo.tif")
    expect_true(file.exists("foo.tif"))
    file.remove("foo.tif")
    file.remove("foo.tif.aux.xml")
})

test_that("Rotate method for SpatRasterImage which converts to ExtImage", {
    suppressWarnings(spi <- SpatRasterImage(rast(img_path)))
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

test_that("scaleImg, SpatRasterImage method, EBI behind the scene", {
    suppressWarnings(spi <- SpatRasterImage(rast(img_path)))
    fct <- 1.5
    spi2 <- scaleImg(spi, fct)
    # OK, right now I'm not sure if it's best to keep the center constant
    # but one can always just translate
    ext_old <- ext(spi)
    ext_new <- ext(spi2)
    expect_equal(bbox_center(ext_old), bbox_center(ext_new))
    old_width <- ext_old["xmax"] - ext_old["xmin"]
    old_height <- ext_old["ymax"] - ext_old["ymin"]
    new_width <- ext_new["xmax"] - ext_new["xmin"]
    new_height <- ext_new["ymax"] - ext_new["ymin"]
    expect_equal(new_width, old_width * fct)
    expect_equal(new_height, old_height * fct)
})

test_that("affineImg, SpatRasterImage method, EBI behind the scene", {
    suppressWarnings(spi <- SpatRasterImage(rast(img_path)))
    # Should rotate clockwise somewhat, squished to be flatter, shear to the right
    M <- matrix(c(0.6, -0.2, 0.2, 0.3), nrow = 2)
    v <- c(0, 300)
    ebi <- affineImg(spi, M, v)
    ext_exp <- ext(spi) |> st_bbox() |> st_as_sfc()
    ext_exp <- st_bbox(ext_exp * t(M) + v) |> as.vector()
    expect_equal(ext_exp, unname(ext(ebi)[c("xmin", "ymin", "xmax", "ymax")]))
    # check content
    expect_equal(imgRaster(ebi)[6,3,2] |> as.numeric(), 191)
    # The resulting image tightly bounds the rotated image
    expect_true(any(imgRaster(ebi)[1,,] > 0))
    expect_true(any(imgRaster(ebi)[,1,] > 0))
    expect_true(any(imgRaster(ebi)[dim(ebi)[1],,] > 0))
    expect_true(any(imgRaster(ebi)[,dim(ebi)[2],] > 0))
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

test_that("Convert SpatRasterImage to ExtImage, RGB", {
    suppressWarnings(spi <- SpatRasterImage(rast(img_path)))
    ebi <- toExtImage(spi)
    expect_s4_class(ebi, "ExtImage")
    expect_equal(ext(spi), ext(ebi))
    expect_equal(imgRaster(spi) |> terra::as.array(),
                 imgRaster(ebi) |> as.array() |> aperm(c(2,1,3)))
})

fp <- tempdir()
vizgen_dir <- VizgenOutput(file_path = file.path(fp, "vizgen_test"))
fn <- file.path(vizgen_dir, "images", "mosaic_Cellbound1_z3.tif")

test_that("Convert SpatRasterImage to ExtImage, grayscale", {
    suppressWarnings(spi <- SpatRasterImage(rast(fn)))
    ebi <- toExtImage(spi)
    expect_equal(ext(spi), ext(ebi))
    m1 <- imgRaster(spi) |> terra::as.array()
    m1 <- m1[,,1]
    expect_equal(m1, imgRaster(ebi) |> as.array() |> t())
})

test_that("Crop SpatRasterImage", {
    suppressWarnings(spi <- SpatRasterImage(rast(fn)))
    bbox_use <- c(xmin = 100, xmax = 200, ymin = 100, ymax = 200)
    spi_sub <- cropImg(spi, bbox_use)
    expect_s4_class(spi_sub, "SpatRasterImage")
    expect_equal(ext(spi_sub), bbox_use)
})

# BioFormatsImage=================
# Often the first time it doesn't work but the second time it does
tryCatch(xml_meta <- read.omexml(xenium_fn) |> read_xml() |> as_list())
xml_meta <- read.omexml(xenium_fn) |> read_xml() |> as_list()
psx <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeX") |> as.numeric()
psy <- attr(xml_meta$OME$Image$Pixels, "PhysicalSizeY") |> as.numeric()
meta <- read.metadata(xenium_fn) |> coreMetadata(1)
sizeX_full <- meta$sizeX
sizeY_full <- meta$sizeY
test_that("BioFormatsImage constructor", {
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
    expect_equal(d, c(X=sizeX_full, Y=sizeY_full, C=1,Z=1,"T"=1))
})

test_that("show for BioFormatsImage", {
    bfi <- BioFormatsImage(xenium_fn)
    expect_output(show(bfi), "X: 1431, Y: 1288, C: 1, Z: 1, T: 1, BioFormatsImage")
})

test_that("Errors when constructing BioFormatsImage", {
    library(RBioFormats)
    # extent
    expect_error(BioFormatsImage(xenium_fn, isFull = FALSE),
                 "Extent must be specified when isFull = FALSE")
    # Invalid bbox
    bbox_use <- c(xmin = 100, xmax = 0, ymin = 0, ymax = 102)
    expect_error(BioFormatsImage(xenium_fn, ext = bbox_use, isFull = FALSE),
                 "Min limit is greater than max limit")
    # isFull
    expect_error(BioFormatsImage(xenium_fn, ext = bbox_use, isFull = NA),
                 "isFull must be either TRUE or FALSE, not NA")
    # origin
    expect_error(BioFormatsImage(xenium_fn, origin = c(1,2,3)),
                 "origin must be a numeric vector of length 2 without NAs.")
    # Transformation
    expect_error(BioFormatsImage(xenium_fn, transformation = list(name = "foo")),
                 "Name of transformation must be one of ")
    expect_error(BioFormatsImage(xenium_fn, transformation = list(M = diag(nrow = 3), v = c(0,0))),
                 "M must be a 2x2 numeric matrix")
})

test_that("Convert BioFormatsImage to ExtImage, full extent", {
    library(RBioFormats)
    ext_expect <- c(xmin = 0, ymin = 0, xmax = sizeX_full*psx, ymax = sizeY_full*psy)
    bfi <- BioFormatsImage(xenium_fn)
    expect_message(toExtImage(bfi, resolution = 10L),
                   "Resolution subscript out of bound")
    ebi <- toExtImage(bfi, resolution = 1L)
    expect_s4_class(ebi, "ExtImage")
    expect_equal(ext(ebi), ext_expect[c("xmin", "xmax", "ymin", "ymax")])
    dim_img <- dim(imgRaster(ebi))
    expect_equal(dim_img, c(sizeX_full, sizeY_full))
})

ext_use <- c(xmin = 1000, xmax = 2000, ymin = 600, ymax = 1600)
test_that("Convert BioFormatsImage to ExtImage, not full extent", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    ebi <- toExtImage(bfi, resolution = 1L)
    expect_s4_class(ebi, "ExtImage")
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

test_that("Convert BioFormatsImage to ExtImage, not resolution 1", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use, isFull = FALSE)
    ebi <- toExtImage(bfi, resolution = 2L)
    expect_s4_class(ebi, "ExtImage")
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
    ext_use <- c(xmin = 0, xmax = 14, ymin = 0, ymax = 23)
    expect_warning(bfi <- BioFormatsImage(img_path),
                   "Physical pixel size absent from image metadata")
    expect_equal(ext(bfi), ext_use)
})

test_that("Ignore resolution in toExtImage when there's only 1 resolution", {
    ext_use <- c(xmin = 0, xmax = 14, ymin = 0, ymax = 23)
    suppressWarnings(bfi <- BioFormatsImage(img_path))
    expect_warning(ebi <- toExtImage(bfi, resolution = 4L),
                   "Physical pixel size absent from image metadata")
    expect_equal(ext(ebi), ext_use)
})

test_that("Convert BioFormatsImage to SpatRasterImage", {
    library(RBioFormats)
    fp <- tempdir()
    fn <- XeniumOutput(file_path = file.path(fp, "xenium_test"))
    fn1 <- file.path(fn, "morphology_mip.ome.tif")
    bfi <- BioFormatsImage(fn1, ext_use, isFull = FALSE)
    expect_message(spi <- toSpatRasterImage(bfi, resolution = 1L), "non OME")
    fn <- file.path(dirname(fn1), "morphology_mip_res1.tiff")
    expect_true(file.exists(fn))
    expect_s4_class(spi, "SpatRasterImage")
    ext_sf_bfi <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf_spi <- st_bbox(ext(spi)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf_bfi, ext_sf_spi, sparse = FALSE))
    expect_true(st_area(ext_sf_spi) / st_area(ext_sf_bfi) < 1.005)
    unlink(fn, recursive = TRUE)
})

test_that("Convert BioFormatsImage to SpatRasterImage not saving geotiff", {
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
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    bfi2 <- transposeImg(bfi, resolution = 1L)
    expect_s4_class(bfi2, "BioFormatsImage")
    expect_equal(transformation(bfi2), list(name="transpose"))
    ext_bfi2 <- ext(bfi2)
    expect_true(ext_bfi2["ymax"] > ext_bfi2["xmax"])
    # Internally the ext is not changed
    expect_equal(bfi@ext, bfi2@ext)
})

test_that("mirror, vertical", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    bfi2 <- mirrorImg(bfi, direction = "vertical")
    expect_s4_class(bfi2, "BioFormatsImage")
    expect_equal(transformation(bfi2), list(name="mirror", direction="vertical"))
    expect_equal(ext(bfi), ext(bfi2))
})

test_that("mirror, horizontal", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    bfi2 <- mirrorImg(bfi, direction = "horizontal")
    expect_s4_class(bfi2, "BioFormatsImage")
    expect_equal(transformation(bfi2), list(name="mirror", direction="horizontal"))
    expect_equal(ext(bfi), ext(bfi2))
})

test_that("rotate BFI", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    bfi2 <- rotateImg(bfi, degrees = 45)
    expect_s4_class(bfi2, "BioFormatsImage")
    expect_equal(transformation(bfi2), list(name="rotate", degrees=45))
    expect_equal(bbox_center(ext(bfi)), bbox_center(ext(bfi2)))
    ext_exp <- st_bbox(ext(bfi)) |> st_as_sfc()
    m <- matrix(c(cos(pi/4), sin(pi/4), -sin(pi/4), cos(pi/4)), ncol = 2)
    cent <- st_centroid(ext_exp) |> st_coordinates()
    ext_exp <- (ext_exp - cent)*m + cent
    expect_equal(st_bbox(ext_exp) |> as.vector(), unname(ext(bfi2)[c("xmin","ymin","xmax","ymax")]))
})

test_that("translateImg, BioFormatsImage method", {
    bfi <- BioFormatsImage(xenium_fn)
    v <- c(135, 246)
    ext_old <- ext(bfi)
    bfi_tr <- translateImg(bfi, v)
    ext_new <- ext(bfi_tr)
    expect_equal(ext_new[c("xmin", "xmax")], ext_old[c("xmin", "xmax")] + v[1])
    expect_equal(ext_new[c("ymin", "ymax")], ext_old[c("ymin", "ymax")] + v[2])
    expect_equal(SpatialFeatureExperiment::origin(bfi_tr), v)
})

test_that("scaleImg, BioFormatsImage method", {
    bfi <- BioFormatsImage(xenium_fn)
    fct <- 1.5
    bfi2 <- scaleImg(bfi, fct)
    # OK, right now I'm not sure if it's best to keep the center constant
    # but one can always just translate
    ext_old <- ext(bfi)
    ext_new <- ext(bfi2)
    expect_equal(bbox_center(ext_old), bbox_center(ext_new))
    old_width <- ext_old["xmax"] - ext_old["xmin"]
    old_height <- ext_old["ymax"] - ext_old["ymin"]
    new_width <- ext_new["xmax"] - ext_new["xmin"]
    new_height <- ext_new["ymax"] - ext_new["ymin"]
    expect_equal(new_width, old_width * fct)
    expect_equal(new_height, old_height * fct)
})

test_that("affineImg, BioFormatsImage method", {
    bfi <- BioFormatsImage(xenium_fn)
    M <- matrix(c(0.6, -0.2, 0.2, 0.3), nrow = 2)
    v <- c(0, 300)
    bfi2 <- affineImg(bfi, M, v)
    expect_equal(transformation(bfi2), list(M=M, v=v))
    ext_exp <- ext(bfi) |> st_bbox() |> st_as_sfc()
    ext_exp <- st_bbox(ext_exp * t(M) + v) |> as.vector()
    expect_equal(ext_exp, unname(ext(bfi2)[c("xmin", "ymin", "xmax", "ymax")]))
})

test_that("Combine transformations for BFI", {
    bfi <- BioFormatsImage(xenium_fn)
    # When they cancel out
    bfi2 <- mirrorImg(bfi, "vertical")
    bfi2 <- mirrorImg(bfi2, "vertical")
    expect_equal(transformation(bfi2), list())
    # Create new matrix
    bfi3 <- mirrorImg(bfi, "vertical")
    bfi3 <- rotateImg(bfi3, degrees = 30)
    expect_named(transformation(bfi3), c("M", "v"))
    expect_equal(bbox_center(ext(bfi)), bbox_center(ext(bfi3)))
    m1 <- matrix(c(1,0,0,-1), ncol = 2)
    m2 <- matrix(c(cos(pi/6), -sin(pi/6), sin(pi/6), cos(pi/6)), ncol = 2)
    expect_equal(transformation(bfi3)$M, m2 %*% m1)
})

test_that("Convert BioFormatsImage to ExtImage after translation", {
    library(RBioFormats)
    # Make sure that the right part of the image is loaded
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- toExtImage(bfi, resolution = 1L)
    v <- c(135, 246)
    bfi_tr <- translateImg(bfi, v)
    ebi_sub <- toExtImage(bfi_tr, resolution = 1L)
    ext_sub <- ext(ebi_sub)
    ext_exp <- ext(ebi)
    ext_exp[c("xmin", "xmax")] <- ext_exp[c("xmin", "xmax")] + v[1]
    ext_exp[c("ymin", "ymax")] <- ext_exp[c("ymin", "ymax")] + v[2]
    expect_equal(ext_sub, ext_exp)
    expect_equal(imgRaster(ebi), imgRaster(ebi_sub))
})

test_that("Convert BFI to EBI after transformation, full extent", {
    bfi <- BioFormatsImage(xenium_fn)
    bfi2 <- transposeImg(bfi)
    ebi <- toExtImage(bfi2, resolution = 1L)
    expect_equal(ext(ebi), ext(bfi2))
    expect_equal(imgRaster(ebi), toExtImage(bfi, resolution = 1L) |> imgRaster() |> EBImage::transpose())
})

test_that("Convert BFI to EBI after transformation, not full extent", {
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    bfi2 <- mirrorImg(bfi, "horizontal")
    ebi <- toExtImage(bfi2, resolution = 1L)
    # ext will be slightly larger due to pixelation and rounding out
    ext_sf1 <- ext(ebi) |> st_bbox() |> st_as_sfc()
    ext_sf2 <- ext(bfi2) |> st_bbox() |> st_as_sfc()
    expect_true(st_covers(ext_sf1, ext_sf2, sparse = FALSE))
    expect_true(st_area(ext_sf1)/st_area(ext_sf2) < 1.005)
    # check content
    expect_equal(imgRaster(ebi), toExtImage(bfi, resolution = 1L) |> imgRaster() |> EBImage::flop())
})

test_that("Crop BioFormatsImage", {
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

# ExtImage====================

test_that("ExtImage constructor", {
    img <- readImage(system.file('images', 'nuclei.tif', package="EBImage"))
    expect_error(ExtImage(img, ext = NULL),
                 "Extent must be specified")
})

test_that("show for ExtImage", {
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi0 <- toExtImage(bfi, resolution = 1L)
    expect_output(show(ebi0), "882 x 588 \\(width x height\\) ExtImage")
})

test_that("transpose, check ext", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi0 <- toExtImage(bfi, resolution = 1L)
    ebi <- transposeImg(ebi0)
    expect_s4_class(ebi, "ExtImage")
    ext_ebi <- ext(ebi)
    expect_true(ext_ebi["ymax"] > ext_ebi["xmax"])
    expect_equal(imgRaster(ebi), ebi0 |> imgRaster() |> EBImage::transpose())
})

test_that("mirror, vertical", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi0 <- toExtImage(bfi, resolution = 1L)
    ebi <- mirrorImg(ebi0, direction = "vertical")
    ext_ebi <- ext(ebi)
    expect_equal(ext_ebi, ext(ebi0))
    expect_equal(imgRaster(ebi), EBImage::flip(imgRaster(ebi0)))
})

test_that("mirror, horizontal", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi0 <- toExtImage(bfi, resolution = 1L)
    ebi <- mirrorImg(ebi0, direction = "horizontal")
    ext_ebi <- ext(ebi)
    expect_equal(ext_ebi, ext(ebi0))
    expect_equal(imgRaster(ebi), EBImage::flop(imgRaster(ebi0)))
})

test_that("translateImg, ExtImage method", {
    library(RBioFormats)
    bfi <- BioFormatsImage(xenium_fn, ext_use2, isFull = FALSE)
    ebi <- toExtImage(bfi, resolution = 1L)
    ext_old <- ext(ebi)
    v <- c(135, 246)
    ebi_tr <- translateImg(ebi, v)
    ext_new <- ext(ebi_tr)
    expect_equal(ext_new[c("xmin", "xmax")], ext_old[c("xmin", "xmax")] + v[1])
    expect_equal(ext_new[c("ymin", "ymax")], ext_old[c("ymin", "ymax")] + v[2])
})

test_that("Crop ExtImage", {
    library(RBioFormats)
    ebi <- BioFormatsImage(xenium_fn) |> toExtImage(resolution = 1L)
    ebi_sub <- cropImg(ebi, ext_use)
    expect_s4_class(ebi_sub, "ExtImage")
    ext_sf1 <- st_bbox(ext_use) |> st_as_sfc()
    ext_sf2 <- st_bbox(ext(ebi_sub)) |> st_as_sfc()
    expect_true(st_covered_by(ext_sf1, ext_sf2, sparse = FALSE))
    expect_true(st_area(ext_sf2) / st_area(ext_sf1) < 1.005)
    # pixel range
    dim_expect <- c(588, 589)
    expect_equal(dim(imgRaster(ebi_sub)), round(dim_expect))
})

# Image setter-------
test_that("Image setter, the image isn't already there", {
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    # Weirdly the first time I get the null pointer error
    sfe <- readXenium(fn)
    img <- getImg(sfe) |> toExtImage(resolution = 1L)
    img <- img > 500
    Img(sfe, image_id = "mask") <- img
    df <- imgData(sfe)
    expect_true("mask" %in% imageIDs(sfe))
    img2 <- getImg(sfe, image_id = "mask")
    expect_equal(img2, img)
    unlink(fn, recursive = TRUE)
})

test_that("Image setter, modify existing image", {
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    # Weirdly the first time I get the null pointer error
    sfe <- readXenium(fn)
    df_old <- imgData(sfe)
    img <- getImg(sfe) |> toExtImage(resolution = 1L)
    # Increase contrast
    img <- img * 1.1
    Img(sfe, image_id = "morphology_focus") <- img
    img2 <- getImg(sfe, image_id = "morphology_focus")
    expect_equal(img2, img)
    df_new <- imgData(sfe)
    expect_equal(nrow(df_new), nrow(df_old))
    unlink(fn, recursive = TRUE)
})

# Final cleanup in case failed test messed with cleanup
fp <- tempdir()
unlink(file.path(fp, "cosmx_test"), recursive = TRUE)
unlink(file.path(fp, "vizgen_test"), recursive = TRUE)
unlink(file.path(fp, "xenium_test"), recursive = TRUE)
