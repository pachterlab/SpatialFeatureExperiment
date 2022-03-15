# Unit test dimGeometry/ies getters and setters
# 5. Get values only for some but not all sample_ids
# 6. Set to NULL
# 8. Set values only for some but not all sample_ids
library(SingleCellExperiment)
library(S4Vectors)
library(sf)
sfe <- readRDS(system.file("testdata/sfe_toy.rds", package = "SpatialFeatureExperiment"))

test_that("Get List of length 0 when dimGeometries are absent", {
  foo <- dimGeometries(sfe, 2)
  expect_true(is(foo, "List"))
  expect_equal(length(foo), 0L)
})

cg_toy <- readRDS(system.file("testdata/cg_toy.rds", package = "SpatialFeatureExperiment"))
cg_toy2 <- readRDS(system.file("testdata/cg_toy2.rds", package = "SpatialFeatureExperiment"))

test_that("colGeometries setter", {
  colGeometries(sfe) <- list(coords = cg_toy, buffered = cg_toy2)
  cgs <- int_colData(sfe)$colGeometries
  expect_true(is(cgs, "DFrame"))
  expect_equal(names(cgs), c("coords", "buffered"))
  expect_true(is(cgs$coords, "sf"))
  expect_true(is(cgs$buffered, "sf"))
})

test_that("colGeometry setter", {
  colGeometry(sfe, "coords") <- cg_toy
  cgs <- int_colData(sfe)$colGeometries
  expect_true(is(cgs, "DFrame"))
  expect_equal(names(cgs), "coords")
  expect_true(is(cgs$coords, "sf"))
})

sfe2 <- sfe
int_colData(sfe2)$colGeometries <- make_zero_col_DFrame(nrow = ncol(sfe2))
int_colData(sfe2)$colGeometries$coords <- cg_toy
int_colData(sfe2)$colGeometries$buffered <- cg_toy2

test_that("colGeometries getter", {
  cgs <- colGeometries(sfe2)
  expect_true(is(cgs, "List"))
  expect_equal(length(cgs), 2L)
  expect_equal(names(cgs), c("coords", "buffered"))
})

test_that("colGeometry getter", {
  cg <- colGeometry(sfe2)
  expect_equal(cg, cg_toy)
  cg2 <- colGeometry(sfe2, "buffered")
  expect_equal(cg2, cg_toy2)
  cg3 <- colGeometry(sfe2, 2L)
  expect_equal(cg3, cg_toy2)
})

test_that("colGeometryNames", {
  nms <- colGeometryNames(sfe2)
  expect_equal(nms, c("coords", "buffered"))
  colGeometryNames(sfe2) <- c("foo", "bar")
  expect_equal(names(int_colData(sfe2)$colGeometries), c("foo", "bar"))
})

# More than one sample_id
sfe3 <- readRDS(system.file("testdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))

test_that("colGeometry setter for one of the two samples (not already present)", {
  # colGeometry not already present
  colGeometry(sfe3, type = "coords", sample_id = "sample01") <- cg_toy[1:3,]
  expect_equal(int_colData(sfe3)$colGeometries$coords[1:3,], cg_toy[1:3,])
  expect_true(all(st_is_empty(int_colData(sfe3)$colGeometries$coords[4:5,])))
})

test_that("colGeometry setter for one of the two samples (already present)", {
  # colGeometry already present
  foo <- st_sf(geometry = st_sfc(lapply(seq_len(ncol(sfe3)),
                                        function(t) st_geometrycollection())),
               sf_column_name = "geometry", row.names = colnames(sfe3))
  int_colData(sfe3)$colGeometries <- make_zero_col_DFrame(nrow = ncol(sfe3))
  int_colData(sfe3)$colGeometries$foo <- foo
  colGeometry(sfe3, "foo", sample_id = "sample01") <- cg_toy[1:3,]
  expect_equal(int_colData(sfe3)$colGeometries$foo[1:3,], cg_toy[1:3,])
  expect_true(all(st_is_empty(int_colData(sfe3)$colGeometries$foo[4:5,])))
})

test_that("colGeometry getter for one of the two samples", {
  int_colData(sfe3)$colGeometries <- make_zero_col_DFrame(nrow = ncol(sfe3))
  int_colData(sfe3)$colGeometries$coords <- cg_toy
  coords_sample02 <- colGeometry(sfe3, "coords", "sample02")
  expect_equal(nrow(coords_sample02), 2L)
  expect_equal(coords_sample02, cg_toy[4:5,])
})
