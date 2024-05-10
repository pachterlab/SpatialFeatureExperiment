outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
samples <- file.path(outdir, paste0("sample0", 1:2))
sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")

test_that("Save SFE with SpatRaster images as RDS", {
    saveRDS(sfe, "foo.rds")
    sfe2 <- readRDS("foo.rds")
    imgs <- imgData(sfe2)$data
    classes <- vapply(imgs, function(x) class(x), FUN.VALUE = character(1))
    expect_true(all(classes == "SpatRasterImage"))
    expect_s4_class(sfe2, "SpatialFeatureExperiment")
    expect_equal(dim(getImg(sfe2)), c(23, 14, 3))
    unlink("foo.rds")
})

test_that("Save SFE without images", {
    imgData(sfe) <- NULL
    saveRDS(sfe, "bar.rds")
    sfe2 <- readRDS("bar.rds")
    expect_equal(sfe, sfe2)
    unlink("bar.rds")
})

test_that("Backward compatibility with old version of SpatRasterImage", {
    fp <- system.file(file.path("extdata", "sfe_old_spi.rds"),
                      package = "SpatialFeatureExperiment")
    sfe <- readRDS(fp)
    img <- getImg(sfe)
    expect_s4_class(img, "SpatRasterImage")
    expect_error(img@image)
    expect_s4_class(img, "SpatRaster")
})
