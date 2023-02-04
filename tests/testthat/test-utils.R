library(SFEData)

sfe <- McKellarMuscleData("small")

test_that("Change sample ID", {
    sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
    expect_equal(sampleIDs(sfe), "sample01")
    expect_equal(names(int_metadata(sfe)$spatialGraphs), "sample01")
    for (n in names(int_metadata(sfe)$annotGeometries)) {
        ag <- int_metadata(sfe)$annotGeometries[[n]]
        expect_equal(unique(ag$sample_id), "sample01")
    }
})
