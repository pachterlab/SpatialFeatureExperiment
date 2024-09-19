sfe <- readRDS(system.file("extdata/sfe_visium.rds", package = "SpatialFeatureExperiment"))
test_that("SFEVersion of first version of SFE object should be NULL", {
    expect_null(SFEVersion(sfe))
})

test_that("Add version in updateObject", {
    curr_version <- packageVersion("SpatialFeatureExperiment")
    expect_message({sfe <- updateObject(sfe, verbose = TRUE)},
                   paste0("Updating it to version ", curr_version))
    expect_equal(SFEVersion(sfe), curr_version)
})

test_that("Update colFeatureData if present", {
    df <- DataFrame(foo = 1, row.names = "sample_id")
    metadata(colData(sfe))$featureData <- df
    sfe <- updateObject(sfe)
    df2 <- colFeatureData(sfe)
    expect_equal(df, df2)
    expect_null(metadata(colData(sfe))$featureData)
})
