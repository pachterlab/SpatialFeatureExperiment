library(SFEData)
library(SpatialExperiment)
library(SingleCellExperiment)
library(terra)
sfe <- McKellarMuscleData("small")
img_path <- system.file(file.path("extdata", "sample01", "outs", "spatial",
                                  "tissue_lowres_image.png"),
                        package = "SpatialFeatureExperiment")
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
})

test_that("transposeImg", {
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    img <- getImg(sfe)
    # SpatRasterImage method
    img_t <- transposeImg(img)
    expect_s4_class(img_t, "SpatRasterImage")
    expect_equal(dim(img@image)[1:2], dim(img_t@image)[2:1])
    # SFE method
    sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    img_t2 <- getImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    expect_equal(dim(img@image)[1:2], dim(img_t2@image)[2:1])
})

test_that("mirrorImg", {
    sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
                  scale_fct = 0.023)
    img <- getImg(sfe)
    # SpatRasterImage method
    img_m <- mirrorImg(img)
    mat1 <- terra::as.array(terra::mean(img@image))[,,1]
    mat2 <- terra::as.array(terra::mean(img_m@image))[,,1]
    mat2_rev <- apply(mat2, 2, rev)
    expect_equal(mat1, mat2_rev)
    # SFE method
    sfe <- mirrorImg(sfe, sample_id = "Vis5A", image_id = "lowres")
    img_m2 <- getImg(sfe)
    mat3 <- terra::as.array(terra::mean(img_m2@image))[,,1]
    mat3_rev <- apply(mat3, 2, rev)
    expect_equal(mat1, mat3_rev)
})

sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres",
              scale_fct = 0.023)
test_that("imgRaster, trivial", {
    img <- imgRaster(getImg(sfe))
    expect_s4_class(img, "SpatRaster")
})

test_that("imgSource, trivial", {
    expect_true(is.na(imgSource(getImg(sfe))))
})
