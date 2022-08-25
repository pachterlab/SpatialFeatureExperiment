library(SingleCellExperiment)
library(S4Vectors)

sfe <- readRDS(system.file("testdata/sfe_toy.rds", package = "SpatialFeatureExperiment"))

test_that("Get List of length 0 when localResults are absent", {
  foo <- localResults(sfe)
  expect_true(is(foo, "List"))
  expect_equal(length(foo), 0L)
})
set.seed(29)
toy_res1 <- matrix(rnorm(10), nrow = 5, ncol = 2,
                   dimnames = list(colnames(sfe), c("meow", "purr")))
toy_res2 <- matrix(rpois(10, lambda = 2), nrow = 5, ncol = 2,
                   dimnames = list(colnames(sfe), c("sassy", "tortitude")))

test_that("localResults setter", {
  localResults(sfe) <- list(foo = toy_res1, bar = toy_res2)
  lrs <- int_colData(sfe)$localResults
  expect_true(is(lrs, "DFrame"))
  expect_equal(names(lrs), c("foo", "bar"))
  expect_equal(lrs$foo, toy_res1)
  expect_equal(lrs$bar, toy_res2)
})

test_that("localResult setter", {
  localResult(sfe, "foo") <- toy_res1
  lrs <- int_colData(sfe)$localResults
  expect_true(is(lrs, "DFrame"))
  expect_equal(names(lrs), "foo")
  expect_equal(lrs$foo, toy_res1)
})

# Should have passed the previous tests
localResults(sfe) <- list(foo = toy_res1, bar = toy_res2)
test_that("localResults getter", {
  lrs <- localResults(sfe)
  expect_true(is(lrs, "List"))
  expect_equal(length(lrs), 2L)
  expect_equal(names(lrs), c("foo", "bar"))
})

test_that("localResults getter", {
  lr <- localResult(sfe)
  expect_equal(lr, toy_res1)
  lr2 <- localResult(sfe, "bar")
  expect_equal(lr2, toy_res2)
  lr3 <- localResult(sfe, 2L)
  expect_equal(lr3, toy_res2)
})

test_that("localResultNames", {
  nms <- localResultNames(sfe)
  expect_equal(nms, c("foo", "bar"))
  localResultNames(sfe) <- c("Moran", "Geary")
  expect_equal(names(int_colData(sfe)$localResults), c("Moran", "Geary"))
})

# More than one sample_id
sfe3 <- readRDS(system.file("testdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))

test_that("localResult setter for one of the two samples (not already present)", {
  localResult(sfe3, type = "foo", sample_id = "sample01",
              withDimnames = FALSE) <- toy_res1[1:3,]
  expect_equal(int_colData(sfe3)$localResults$foo[1:3,], toy_res1[1:3,])
  expect_true(all(is.na(int_colData(sfe3)$localResults$foo[4:5,])))
})

test_that("localResult setter for one of the two samples (already present)", {
  localResults(sfe3) <- list(bar = toy_res2)
  colnames(toy_res1) <- colnames(toy_res2)
  localResult(sfe3, "bar", sample_id = "sample01") <- toy_res1[1:3,]
  expect_equal(int_colData(sfe3)$localResults$bar[1:3,], toy_res1[1:3,],
               ignore_attr = "row.names")
  expect_equal(int_colData(sfe3)$localResults$bar[4:5,], toy_res2[4:5,],
               ignore_attr = "row.names")
})

test_that("localResult getter for one of the two samples", {
  localResults(sfe3) <- list(bar = toy_res2)
  coords_sample02 <- localResult(sfe3, "bar", "sample02")
  expect_equal(nrow(coords_sample02), 2L)
  expect_equal(coords_sample02, toy_res2[4:5,], ignore_attr = "row.names")
})
