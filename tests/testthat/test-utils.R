library(SFEData)
library(sf)

sfe <- McKellarMuscleData("small")
rg <- matrix(rnorm(2*nrow(sfe)), ncol = 2)
colnames(rg) <- c("x", "y")
rg <- as.data.frame(rg)
rg <- st_as_sf(rg, coords = c("x", "y"), crs = NA,
               row.names = rownames(sfe))
rowGeometry(sfe, "foo_Vis5A") <- rg

test_that("Change sample ID", {
    sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
    expect_equal(sampleIDs(sfe), "sample01")
    expect_equal(names(int_metadata(sfe)$spatialGraphs), "sample01")
    for (n in names(int_metadata(sfe)$annotGeometries)) {
        ag <- int_metadata(sfe)$annotGeometries[[n]]
        expect_equal(unique(ag$sample_id), "sample01")
    }
    expect_equal(rowGeometryNames(sfe), "foo_sample01")
})
