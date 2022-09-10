library(sf)
library(SingleCellExperiment)

sfe <- readRDS(system.file("extdata/sfe_toy.rds",
                           package = "SpatialFeatureExperiment"))
ag <- readRDS(system.file("extdata/ag.rds",
                          package = "SpatialFeatureExperiment"))
test_that("annotGeometries setter", {
  annotGeometries(sfe) <- list(hull = ag)
  expect_true(is.list(int_metadata(sfe)$annotGeometries))
  expect_equal(names(int_metadata(sfe)$annotGeometries), "hull")
  expect_equal(int_metadata(sfe)$annotGeometries[["hull"]], ag)
  # Error if can't be converted to sf
  expect_error(annotGeometries(sfe) <- list(foo = iris))
  # Error if sample_id is absent from the sfe object
  foo <- ag
  foo$sample_id <- "bar"
  expect_error(annotGeometries(sfe) <- list(foo = foo))
})

test_that("annotGeometries getter", {
  int_metadata(sfe)$annotGeometries <- list(hull = ag)
  out <- annotGeometries(sfe)
  expect_true(is.list(out))
  expect_equal(names(out), "hull")
  expect_equal(out[["hull"]], ag)
})

test_that("annotGeometry setter (one sample_id)", {
  annotGeometry(sfe, "hull") <- ag
  expect_true(is.list(int_metadata(sfe)$annotGeometries))
  expect_equal(names(int_metadata(sfe)$annotGeometries), "hull")
  expect_equal(int_metadata(sfe)$annotGeometries[["hull"]], ag)
})

sfe3 <- sfe
buffered <- st_buffer(ag, dist = 0.1)
int_metadata(sfe3)$annotGeometries <- list(hull = ag, buffered = buffered)

test_that("annotGeometry getter (one sample_id)", {
  out <- annotGeometry(sfe3, type = "buffered")
  expect_true(is(out, "sf"))
  expect_equal(out, buffered)
  out <- annotGeometry(sfe3)
  expect_equal(out, ag)
})

test_that("annotGeometryNames getter", {
  expect_equal(annotGeometryNames(sfe3), c("hull", "buffered"))
})

test_that("annotGeometryNames setter", {
  annotGeometryNames(sfe3) <- c("foo", "bar")
  expect_equal(names(int_metadata(sfe3)$annotGeometries), c("foo", "bar"))
})

# More than one sample_id
sfe2 <- readRDS(system.file("extdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))
ag2 <- readRDS(system.file("extdata/ag_samples.rds",
                           package = "SpatialFeatureExperiment"))
test_that("annotGeometry getter for one out of two sample_id", {
  int_metadata(sfe2)$annotGeometries <- list(annot = ag2)
  out <- annotGeometry(sfe2, "annot", sample_id = "sample01")
  expect_equal(out, ag2[1,])
})

test_that("annotGeometry setter for one sample_id when already present", {
  # when annotGeometry of the given name already exists
  foo <- st_sf(geometry = st_sfc(lapply(1:2, function(t) st_geometrycollection())),
               sample_id = c("sample01", "sample02"),
               sf_column_name = "geometry")
  int_metadata(sfe2)$annotGeometries$foo <- foo
  annotGeometry(sfe2, "foo", sample_id = "sample01") <- ag2[1,]
  bar <- int_metadata(sfe2)$annotGeometries$foo
  expect_equal(bar[bar$sample_id == "sample01",], ag2[1,])
})

# How about adding annotation for one sample when another sample is already present?
