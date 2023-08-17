outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
samples <- file.path(outdir, paste0("sample0", 1:2))
sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")

test_that("Save SFE with SpatRaster images as RDS", {
    saveRDS(sfe, "foo.rds")
    sfe2 <- readRDS("foo.rds")
    imgs <- imgData(sfe2)$data
    classes <- vapply(imgs, function(x) class(x@image), FUN.VALUE = character(1))
    expect_true(all(classes == "PackedSpatRaster"))
    expect_s4_class(sfe2, "SpatialFeatureExperiment")
    unlink("foo.rds")
})

test_that("Save SFE without images", {
    imgData(sfe) <- NULL
    saveRDS(sfe, "bar.rds")
    sfe2 <- readRDS("bar.rds")
    expect_equal(sfe, sfe2)
    unlink("bar.rds")
})
