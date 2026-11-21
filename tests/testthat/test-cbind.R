sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
    package = "SpatialFeatureExperiment"
))
g_visium <- readRDS(system.file("extdata/colgraph_visium.rds",
    package = "SpatialFeatureExperiment"
))
g_visium2 <- readRDS(system.file("extdata/colgraph_visium2.rds",
    package = "SpatialFeatureExperiment"
))
colGraph(sfe_visium, "visium1", "sample01") <- g_visium
colGraph(sfe_visium, "visium2", "sample02") <- g_visium2
ag_samples <- readRDS(system.file("extdata/ag_samples.rds",
    package = "SpatialFeatureExperiment"
))
annotGeometry(sfe_visium, "foo", sample_id = "all") <- ag_samples
sfe_visium1 <- sfe_visium[, colData(sfe_visium)$sample_id == "sample01"]
sfe_visium2 <- sfe_visium[, colData(sfe_visium)$sample_id == "sample02"]

sfe_concat <- cbind(sfe_visium1, sfe_visium2)

test_that("cbind properly deals with annotGeometries and spatialGraphs", {
    ag_orig <- annotGeometries(sfe_visium)
    ag_concat <- annotGeometries(sfe_concat)
    ag_concat <- ag_concat[names(ag_orig)]
    expect_equal(ag_orig, ag_concat)
    expect_equal(spatialGraphs(sfe_visium), spatialGraphs(sfe_concat))
})

cg_toy <- readRDS(system.file("extdata/cg_toy.rds",
                              package = "SpatialFeatureExperiment"
))
cg_toy2 <- readRDS(system.file("extdata/cg_toy2.rds",
                               package = "SpatialFeatureExperiment"
))
library(DuckDBSpatial)
library(sf)
test_that("cbind, with DuckDBSpatial colGeometries", {
    fn1 <- tempfile(fileext = ".parquet")
    fn2 <- tempfile(fileext = ".parquet")
    writeGeoParquet(cg_toy, fn1)
    writeGeoParquet(cg_toy2[c(1:5,1:3),], fn2)
    cg1 <- readGeoParquet(fn1)
    cg2 <- readGeoParquet(fn2)
    colGeometry(sfe_visium1, "foo") <- cg1
    colGeometry(sfe_visium2, "foo") <- cg2
    sfe_c <- cbind(sfe_visium1, sfe_visium2)
    cg <- colGeometry(sfe_c, "foo", sample_id = "all")
    expect_s4_class(cg, "DuckDBDataFrame")
    expected <- c(cg_toy$geometry, cg_toy2$geometry[c(1:5,1:3)])
    actual <- st_as_sf(cg, crs = NA)
    expect_equal(expected, actual$geometry)
    file.remove(fn1, fn2)
})

# With rowGeometries
rg_toy1 <- cg_toy[1:2,]
rownames(rg_toy1) <- rownames(sfe_visium)
rg_toy2 <- cg_toy[3:4,]
rownames(rg_toy2) <- rownames(sfe_visium)
rg_toy3 <- cg_toy[4:5,]
rownames(rg_toy3) <- rownames(sfe_visium)

test_that("Deal with duplicate sample_id", {
    sfe_visium2 <- changeSampleIDs(sfe_visium2, c(sample02 = "sample01"))
    # I need to check what's so slow in cbind but not the most urgent
    expect_message(sfe2 <- cbind(sfe_visium1, sfe_visium2),
                   "'sample_id's are duplicated across")
    expect_setequal(sampleIDs(sfe2), c("sample01", "sample01_1"))
})

test_that("rowGeometry, both x and y has 1 sample, rowGeometry names don't have samples", {
    txSpots(sfe_visium1) <- rg_toy1
    txSpots(sfe_visium2) <- rg_toy2
    sfe2 <- cbind(sfe_visium1, sfe_visium2)
    expect_equal(rowGeometryNames(sfe2), c("txSpots_sample01", "txSpots_sample02"))
    expect_equal(txSpots(sfe2, "sample01"), rg_toy1)
    expect_equal(txSpots(sfe2, "sample02"), rg_toy2)
})

test_that("rowGeometry, multiple samples", {
    txSpots(sfe_visium, "sample01") <- rg_toy1
    txSpots(sfe_visium, "sample02") <- rg_toy2
    txSpots(sfe_visium, "all") <- rg_toy3
    sfe_visiumb <- sfe_visium
    sfe_visiumb <- changeSampleIDs(sfe_visiumb, c(sample01 = "foo", sample02 = "bar"))
    sfe2 <- cbind(sfe_visium, sfe_visiumb)
    expect_equal(rowGeometryNames(sfe2),
                 c("txSpots_sample01", "txSpots_sample02", "txSpots", "txSpots_foo",
                   "txSpots_bar", "txSpots_1"))
    expect_equal(txSpots(sfe2, "sample01"), rg_toy1)
    expect_equal(txSpots(sfe2, "all"), rg_toy3)
    expect_equal(txSpots(sfe2, "bar"), rg_toy2)
})

test_that("What if some but not all of the objects don't have rowGeometries", {
    txSpots(sfe_visium1) <- rg_toy1
    sfe2 <- cbind(sfe_visium1, sfe_visium2)
    expect_equal(rowGeometryNames(sfe2), "txSpots_sample01")
    expect_equal(txSpots(sfe2, "sample01"), rg_toy1)
})

test_that("cbind when only one or no SFE object is specified", {
    expect_null(SpatialFeatureExperiment::cbind())
    sfe2 <- cbind(sfe_visium1)
    expect_equal(sfe2, sfe_visium1)
})

test_that("cbind more than 2 SFE objects", {
    sfe_visium3 <- changeSampleIDs(sfe_visium1, c("sample01" = "sample03"))
    colnames(sfe_visium3) <- paste0(colnames(sfe_visium3), "-3")
    sfe_c <- cbind(sfe_visium1, sfe_visium2, sfe_visium3)

    ag_orig <- annotGeometries(sfe_visium)
    ag_concat <- annotGeometries(sfe_c)
    ag_concat <- ag_concat[names(ag_orig)]
    expect_equal(ag_orig$foo, ag_concat$foo[ag_concat$foo$sample_id!="sample03",])
    expect_equal(ag_orig$foo[ag_orig$foo$sample_id == "sample01","geometry"],
                 ag_concat$foo[ag_concat$foo$sample_id == "sample03","geometry"])
    expect_equal(spatialGraphs(sfe_visium), spatialGraphs(sfe_c)[c("sample01", "sample02")])
    expect_equal(spatialGraphs(sfe_visium)[["sample01"]], spatialGraphs(sfe_c)[["sample03"]])
})
