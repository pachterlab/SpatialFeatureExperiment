library(SFEData)
library(sf)
library(scater)
#library(Voyager)
devtools::load_all("~/Voyager/")
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
