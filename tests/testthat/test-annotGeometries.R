library(sf)
library(SingleCellExperiment)
library(DuckDBSpatial)

sfe <- readRDS(system.file("extdata/sfe_toy.rds",
    package = "SpatialFeatureExperiment"
))
ag <- readRDS(system.file("extdata/ag.rds",
    package = "SpatialFeatureExperiment"
))
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

test_that("annotGeometries setter, DuckDBDataFrame", {
    fn <- tempfile(fileext = ".parquet")
    writeGeoParquet(ag, fn)
    ag2 <- readGeoParquet(fn)
    annotGeometries(sfe) <- list(hull = ag2)
    expect_s4_class(int_metadata(sfe)$annotGeometries[["hull"]], "DuckDBDataFrame")
    file.remove(fn)
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

test_that("annotGeometry setter (DuckDB, 1 sample)", {
    int_metadata(sfe)$annotGeometries <- NULL
    fn <- tempfile(fileext = ".parquet")
    writeGeoParquet(ag, fn)
    ag2 <- readGeoParquet(fn)
    annotGeometry(sfe, "hull") <- ag2
    expect_s4_class(int_metadata(sfe)$annotGeometries[["hull"]], "DuckDBDataFrame")
    file.remove(fn)
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
    package = "SpatialFeatureExperiment"
))
ag2 <- readRDS(system.file("extdata/ag_samples.rds",
    package = "SpatialFeatureExperiment"
))
test_that("annotGeometry getter for one out of two sample_id", {
    int_metadata(sfe2)$annotGeometries <- list(annot = ag2)
    out <- annotGeometry(sfe2, "annot", sample_id = "sample01")
    expect_equal(out, ag2[1, ])
})

test_that("Error message when the geometry of interest is absent", {
    expect_error(annotGeometry(sfe3, "meow"), "annotGeometry meow is absent")
})

test_that("annotGeometry setter for one sample_id when already present", {
    # when annotGeometry of the given name already exists
    foo <- st_as_sf(data.frame(x=c(0,1), y=c(0,1), sample_id = c("sample01", "sample02")),
                    coords = c("x", "y"))
    int_metadata(sfe2)$annotGeometries$foo <- foo
    annotGeometry(sfe2, "foo", sample_id = "sample01") <- ag2[1, ]
    bar <- int_metadata(sfe2)$annotGeometries$foo
    expect_equal(bar[bar$sample_id == "sample01", ], ag2[1, ])
})
# Tests about DuckDB
test_that("Error message when class of value doesn't match that of existing", {
    fn <- tempfile(fileext = ".parquet")
    writeGeoParquet(ag2[1,], fn)
    ag2_duck <- readGeoParquet(fn)
    int_metadata(sfe2)$annotGeometries$foo <- ag2[2,]
    expect_error(annotGeometry(sfe2, "foo", sample_id = "sample01") <- ag2_duck,
                 "value must have class matching existing entry")
    file.remove(fn)
})

test_that("set for one sample out of two, DuckDB", {
    # more general case where cols don't match
    fn1 <- tempfile(fileext = ".parquet")
    fn2 <- tempfile(fileext = ".parquet")
    writeGeoParquet(ag2[1,], fn1)
    bar <- ag2[2,]
    bar$foobar <- "blah"
    writeGeoParquet(bar, fn2)
    ag_duck1 <- readGeoParquet(fn1)
    ag_duck2 <- readGeoParquet(fn2)
    int_metadata(sfe2)$annotGeometries$foo <- ag_duck1
    annotGeometry(sfe2, "foo", "sample02") <- ag_duck2
    actual_duck <- int_metadata(sfe2)$annotGeometries$foo
    expect_equal(as.vector(actual_duck$foobar) |> unname(), c(NA, "blah"))
    actual_geom <- st_as_binary(actual_duck)
    actual_geom <- dplyr::collect(actual_geom@conn |> dplyr::select(geometry)) |> st_as_sf()
    eq <- st_equals(actual_geom$geometry, ag2$geometry) |> unclass() |> unlist()
    expect_equal(eq, seq_len(nrow(ag2)))
    file.remove(fn1, fn2)
})
