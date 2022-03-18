# Unit test spatialGraphs getters and setters
# 8. 1-5, for getters
# 9. Set to NULL
library(SingleCellExperiment)
sfe2 <- readRDS(system.file("testdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))
cgr1 <- readRDS(system.file("testdata/colgraph1.rds",
                            package = "SpatialFeatureExperiment"))
cgr2 <- readRDS(system.file("testdata/colgraph2.rds",
                            package = "SpatialFeatureExperiment"))
agr1 <- readRDS(system.file("testdata/annotgraph1.rds",
                            package = "SpatialFeatureExperiment"))
agr2 <- readRDS(system.file("testdata/annotgraph2.rds",
                            package = "SpatialFeatureExperiment"))

test_that("Set all applicable sample_ids and margins", {
  spatialGraphs(sfe2) <- list(sample01 = list(col = list(foo = cgr1),
                                              annot = list(bar = agr1)),
                              sample02 = list(col = list(bar = cgr2),
                                              annot = list(foo = agr2)))
  df <- int_metadata(sfe2)$spatialGraphs
  expect_true(is(df, "DFrame"))
  expect_equal(names(df), c("sample01", "sample02"))
  expect_equal(rownames(df), c("row", "col", "annot"))
  expect_true(is.null(df[1,1][[1]]))
  expect_equal(df[2,1][[1]]$foo, cgr1)
  expect_equal(df[2,2][[1]]$bar, cgr2)
  expect_equal(df[3,1][[1]]$bar, agr1)
  expect_equal(df[3,2][[1]]$foo, agr2)
})

test_that("Set all applicable sample_ids in one margin", {
  # spatialGraphs absent
  spatialGraphs(sfe2, MARGIN = 2L) <- list(sample01 = list(foo = cgr1, bar = cgr1),
                                           sample02 = list(baz = cgr2))
  df <- int_metadata(sfe2)$spatialGraphs
  expect_true(is(df, "DFrame"))
  expect_equal(names(df), c("sample01", "sample02"))
  expect_equal(rownames(df), c("row", "col", "annot"))
  expect_true(is.null(df[1,1][[1]]))
  expect_equal(df[2,1][[1]]$foo, cgr1)
  expect_equal(df[2,1][[1]]$bar, cgr1)
  expect_equal(df[2,2][[1]]$baz, cgr2)
  # spatialGraphs present
  spatialGraphs(sfe2, MARGIN = 3L) <- list(sample02 = list(foo = agr2))
  df <- int_metadata(sfe2)$spatialGraphs
  expect_equal(df[3,"sample02"][[1]]$foo, agr2)
})

test_that("Set all applicable margins in one sample_id", {
  spatialGraphs(sfe2, sample_id = "sample01") <- list(col = list(foo = cgr1),
                                                      annot = list(bar = agr1,
                                                                   baz = agr1))
  df <- int_metadata(sfe2)$spatialGraphs
  expect_true(is(df, "DFrame"))
  expect_equal(names(df), c("sample01", "sample02"))
  expect_equal(rownames(df), c("row", "col", "annot"))
  expect_true(is.null(df[1,1][[1]]))
  expect_equal(df[2,1][[1]]$foo, cgr1)
  expect_equal(df[3,1][[1]]$bar, agr1)
  expect_equal(df[3,1][[1]]$baz, agr1)
})

test_that("Set all items in one margin and one sample_id", {
  spatialGraphs(sfe2, MARGIN = 2L, sample_id = "sample02") <- list(foo = cgr2, bar = cgr2)
  df <- int_metadata(sfe2)$spatialGraphs
  expect_true(is(df, "DFrame"))
  expect_equal(names(df), c("sample01", "sample02"))
  expect_equal(rownames(df), c("row", "col", "annot"))
  expect_true(is.null(df[2,1][[1]]))
  expect_equal(df[2,2][[1]]$foo, cgr2)
  expect_equal(df[2,2][[1]]$bar, cgr2)
})

test_that("Set one item in one margin and one sample_id", {
  spatialGraph(sfe2, "foo", 3L, "sample01") <- agr1
  df <- int_metadata(sfe2)$spatialGraphs
  expect_true(is(df, "DFrame"))
  expect_equal(names(df), c("sample01", "sample02"))
  expect_equal(rownames(df), c("row", "col", "annot"))
  expect_true(is.null(df[2,1][[1]]))
  expect_equal(df[3,1][[1]]$foo, agr1)
  # Error when wrong length
  expect_error(spatialGraph(sfe2, "foo", 2L, "sample01") <- cgr2,
               "number of samples")
  # Error when not listw
  expect_error(spatialGraph(sfe2, "foo", 3L, "sample01") <- "foobar",
               "listw")
})
