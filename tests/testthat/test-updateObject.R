sfe <- readRDS(system.file("extdata/sfe_visium.rds", package = "SpatialFeatureExperiment"))
# All hell broke loose after Herve ran updateObject::updatePackageObject()
# and I didn't read the commit message after pulling from upstream devel
test_that("Add version in updateObject", {
    curr_version <- packageVersion("SpatialFeatureExperiment")
    int_metadata(sfe)$SFE_version <- NULL
    expect_message({sfe <- updateObject(sfe, verbose = TRUE)},
                   paste0("Updating it to version ", curr_version))
    expect_equal(SFEVersion(sfe), curr_version)
})

test_that("Update colFeatureData if present", {
    int_metadata(sfe)$SFE_version <- NULL
    df <- DataFrame(foo = 1, row.names = "sample_id")
    metadata(colData(sfe))$featureData <- df
    sfe <- updateObject(sfe)
    df2 <- colFeatureData(sfe)
    expect_equal(df, df2)
    expect_null(metadata(colData(sfe))$featureData)
})
