library(SFEData)
library(sf)
library(scater)
library(Voyager)

sfe <- McKellarMuscleData("small")
rg <- matrix(rnorm(2*nrow(sfe)), ncol = 2)
colnames(rg) <- c("x", "y")
rg <- as.data.frame(rg)
rg <- st_as_sf(rg, coords = c("x", "y"), crs = NA,
               row.names = rownames(sfe))
rowGeometry(sfe, "foo_Vis5A") <- rg

sfe <- sfe[,sfe$in_tissue]
colGraph(sfe, "visium") <- findVisiumGraph(sfe, style = "S")
sfe <- logNormCounts(sfe)
sfe <- runPCA(sfe, ncomponents = 10)
annotGraph(sfe, "myofiber_tri2nb") <-
    findSpatialNeighbors(sfe, type = "myofiber_simplified", MARGIN = 3L,
                         method = "tri2nb", dist_type = "idw",
                         zero.policy = TRUE
    )

test_that("Change sample ID", {
    sfe <- runMoransI(sfe, rownames(sfe)[1])
    sfe <- colDataMoransI(sfe, "nCounts")
    sfe <- annotGeometryMoransI(sfe, features = "area",
                                annotGraphName = "myofiber_tri2nb",
                                annotGeometryName = "myofiber_simplified",
                                zero.policy = TRUE
    )
    sfe <- reducedDimMoransI(sfe, "PCA", components = 1:3)
    sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
    expect_equal(sampleIDs(sfe), "sample01")
    expect_equal(names(int_metadata(sfe)$spatialGraphs), "sample01")
    for (n in names(int_metadata(sfe)$annotGeometries)) {
        ag <- int_metadata(sfe)$annotGeometries[[n]]
        expect_equal(unique(ag$sample_id), "sample01")
    }
    expect_equal(rowGeometryNames(sfe), "foo_sample01")

    name_expect <- c("moran_sample01", "K_sample01")
    cfd <- colFeatureData(sfe)
    expect_equal(names(cfd), name_expect)
    afd <- geometryFeatureData(sfe, "myofiber_simplified", 3L)
    expect_equal(names(afd), name_expect)
    rfd <- reducedDimFeatureData(sfe, "PCA")
    expect_equal(names(rfd), name_expect)
})

test_that("bbox center", {
    bbox_use <- c(xmin = 10, xmax = 20, ymin = 5, ymax = 15)
    expect_equal(bbox_center(bbox_use), c(15, 10))
})
fp <- tempdir()
dir_use <- XeniumOutput("v1", file_path = file.path(fp, "xenium_test"))
test_that("Get pixel size", {
    library(RBioFormats)
    fn <- file.path(dir_use, "morphology_focus.ome.tif")
    # Top resolution
    try(pss <- getPixelSize(fn))
    pss <- getPixelSize(fn)
    expect_equal(pss, c(1.7, 1.7))
    pss2 <- getPixelSize(fn, 2L)
    expect_equal(pss2, c(5.1, 5.1))
})

test_that("Aggregate bboxes", {
    bboxes <- list(c(xmin = 5, xmax = 10, ymin = 2, ymax = 20),
                   c(xmin = 8, xmax = 18, ymin = 0, ymax = 15))
    bbox_all <- aggBboxes(bboxes)
    expect_equal(bbox_all, c(xmin = 5, ymin = 0, xmax = 18, ymax = 20))
})

test_that("Image IDs", {
    sfe <- readXenium(dir_use)
    expect_setequal(imageIDs(sfe), c("morphology_focus", "morphology_mip"))
})

unlink(dir_use, recursive = TRUE)
