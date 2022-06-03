# Unit test validity functions
library(sf)
library(SingleCellExperiment)

sfe <- readRDS(system.file("testdata/sfe_toy.rds", package = "SpatialFeatureExperiment"))
cg_toy <- readRDS(system.file("testdata/cg_toy.rds", package = "SpatialFeatureExperiment"))
cg_toy2 <- readRDS(system.file("testdata/cg_toy2.rds", package = "SpatialFeatureExperiment"))
ag <- readRDS(system.file("testdata/ag.rds", package = "SpatialFeatureExperiment"))
cgr1 <- readRDS(system.file("testdata/colgraph1.rds",
                            package = "SpatialFeatureExperiment"))

test_that("Everything in *Geometries must be sf objects", {
  int_colData(sfe)$colGeometries <- make_zero_col_DFrame(nrow = ncol(sfe))
  int_colData(sfe)$colGeometries$coords <- cg_toy
  expect_true(validObject(sfe))
  int_colData(sfe)$colGeometries$coords <- st_drop_geometry(cg_toy)
  expect_error(validObject(sfe), "rather\\s+than\\s+sf")
})

test_that("All sample_ids in annotGeometries must be present in colData", {
  int_metadata(sfe)$annotGeometries$foo <- ag
  expect_true(validObject(sfe))
  ag$sample_id <- "foo"
  int_metadata(sfe)$annotGeometries$foo <- ag
  expect_error(validObject(sfe), "absent from colData")
  ag$sample_id <- NULL
  int_metadata(sfe)$annotGeometries$foo <- ag
  expect_error(validObject(sfe), "does not have column sample_id")
})

test_that("All graphs must be listw objects", {
  int_metadata(sfe)$spatialGraphs <- SpatialFeatureExperiment:::.initialize_spatialGraphs(sfe)
  int_metadata(sfe)$spatialGraphs[3,1] <- I(list(list(foo = "bar")))
  expect_error(validObject(sfe), "must have class listw")
})

test_that("All row and col graphs must have the right number of nodes", {
  int_metadata(sfe)$spatialGraphs <- SpatialFeatureExperiment:::.initialize_spatialGraphs(sfe)
  int_metadata(sfe)$spatialGraphs[2,1] <- I(list(list(foo = cgr1)))
  expect_error(validObject(sfe), "do not have the right length")
})

test_that("The spatialGraphs field must have the right structure", {
  int_metadata(sfe)$spatialGraphs <- list(foo = "bar")
  expect_error(validObject(sfe), "must be a DataFrame")
  df <- SpatialFeatureExperiment:::.initialize_spatialGraphs(sfe)
  names(df) <- "foo"
  int_metadata(sfe)$spatialGraphs <- df
  expect_error(validObject(sfe), "but not colData")
  names(df) <- "sample01"
  rownames(df) <- c("foo", "bar", "baz")
  int_metadata(sfe)$spatialGraphs <- df
  expect_error(validObject(sfe), "Row names")
})
