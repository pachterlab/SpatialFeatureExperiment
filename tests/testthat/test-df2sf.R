library(sf)
pts_df <- readRDS(system.file("extdata/pts_df.rds",
    package = "SpatialFeatureExperiment"
))
pts_sf <- readRDS(system.file("extdata/pts_sf.rds",
    package = "SpatialFeatureExperiment"
))
test_that("df2sf works properly for points", {
    sf_use <- df2sf(pts_df, geometryType = "POINT")
    expect_equal(sf_use, pts_sf, ignore_attr = TRUE)
})

test_that("Points but with spotDiameter", {
    sf_use <- df2sf(pts_df, geometryType = "POINT", spotDiameter = 0.1)
    pts_sf_dia <- readRDS(system.file("extdata/pts_sf_dia.rds",
        package = "SpatialFeatureExperiment"
    ))
    expect_equal(sf_use, pts_sf_dia, ignore_attr = TRUE)
})

test_that("sample_id when I need to split the data frame", {
    multipts_df <- readRDS(system.file("extdata/multipts_df.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- df2sf(multipts_df, geometryType = "MULTIPOINT")
    multipts_sf <- readRDS(system.file("extdata/multipts_sf.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- sf_use[, names(multipts_sf)]
    expect_equal(sf_use, multipts_sf, ignore_attr = TRUE)

    multipts_df_wrong_sample <- readRDS(system.file("extdata/pts_sf_dia.rds",
        package = "SpatialFeatureExperiment"
    ))
    expect_error(df2sf(multipts_df_wrong_sample, geometryType = "MULTIPOINT"))
})

test_that("Linestring", {
    ls_df <- readRDS(system.file("extdata/ls_df.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- df2sf(ls_df, geometryType = "LINESTRING")
    ls_sf <- readRDS(system.file("extdata/ls_sf.rds",
        package = "SpatialFeatureExperiment"
    ))
    expect_equal(sf_use, ls_sf, ignore_attr = TRUE)
})

test_that("Not enough vertices for the specified geometry", {
    ls_df_singleton <- readRDS(system.file("extdata/ls_df_singleton.rds",
        package = "SpatialFeatureExperiment"
    ))
    expect_warning(
        sf_use <- df2sf(ls_df_singleton, geometryType = "LINESTRING"),
        "Removed"
    )
    ls_sf_singleton <- readRDS(system.file("extdata/ls_sf_singleton.rds",
        package = "SpatialFeatureExperiment"
    ))
    expect_equal(sf_use, ls_sf_singleton, ignore_attr = TRUE)
})

test_that("Multilinestring", {
    multils_df <- readRDS(system.file("extdata/multils_df.rds",
        package = "SpatialFeatureExperiment"
    ))
    multils_sf <- readRDS(system.file("extdata/multils_sf.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- df2sf(multils_df, geometryType = "MULTILINESTRING")
    expect_equal(multils_sf, sf_use, ignore_attr = TRUE)
})

test_that("Polygons", {
    pol_df <- readRDS(system.file("extdata/pol_df.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- df2sf(pol_df,
        geometryType = "POLYGON",
        spatialCoordsNames = c("V1", "V2")
    )
    pol_sf <- readRDS(system.file("extdata/pol_sf.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- sf_use[, names(pol_sf)]
    expect_equal(pol_sf, sf_use, ignore_attr = TRUE)
})

test_that("De facto points get POINT geometry", {
    sf_use <- df2sf(pts_df, geometryType = "POLYGON")
    expect_true(all(st_is(sf_use, "POINT")))
    expect_equal(sf_use, pts_sf, ignore_attr = TRUE)
})

test_that("Multipolygons", {
    mpol_df <- readRDS(system.file("extdata/mpol_df.rds",
        package = "SpatialFeatureExperiment"
    ))
    mpol_sf <- readRDS(system.file("extdata/mpol_sf.rds",
        package = "SpatialFeatureExperiment"
    ))
    sf_use <- df2sf(mpol_df,
        geometryType = "MULTIPOLYGON",
        spatialCoordsNames = c("V1", "V2")
    )
    sf_use <- sf_use[, names(mpol_sf)]
    expect_equal(sf_use, mpol_sf, ignore_attr = TRUE)
})
