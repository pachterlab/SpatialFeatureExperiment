# Unit test spatialGraphs getters and setters
library(SingleCellExperiment)
sfe2 <- readRDS(system.file("extdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))
cgr1 <- readRDS(system.file("extdata/colgraph1.rds",
                            package = "SpatialFeatureExperiment"))
cgr2 <- readRDS(system.file("extdata/colgraph2.rds",
                            package = "SpatialFeatureExperiment"))
agr1 <- readRDS(system.file("extdata/annotgraph1.rds",
                            package = "SpatialFeatureExperiment"))
agr2 <- readRDS(system.file("extdata/annotgraph2.rds",
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

test_that("Set graph of the same name for multiple sample_ids", {
  value <- list(sample01 = cgr1, sample02 = cgr2)
  spatialGraphs(sfe2, MARGIN = 2L, sample_id = "all", name = "foo") <- value
  expect_equal(spatialGraphNames(sfe2, 2, "sample01"), "foo")
  expect_equal(spatialGraphNames(sfe2, 2, "sample02"), "foo")
  expect_equal(spatialGraph(sfe2, "foo", 2, "sample01"), cgr1)
  expect_equal(spatialGraph(sfe2, "foo", 2, "sample02"), cgr2)
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

# If it passes the above tests, then the setters should work. Then I'll use the
# setter rather than to use int_metadata to unit test the getters.
spatialGraphs(sfe2) <- list(sample01 = list(col = list(foo = cgr1),
                                            annot = list(bar = agr1)),
                            sample02 = list(col = list(bar = cgr2),
                                            annot = list(foo = agr2)))
mar_names <- c("row", "col", "annot")
test_that("Get all applicable sample_ids and margins", {
  out_all <- spatialGraphs(sfe2)
  expect_equal(names(out_all), c("sample01", "sample02"))
  expect_equal(names(out_all[[1]]), mar_names)
  expect_equal(names(out_all[[2]]), mar_names)
  expect_true(is.null(out_all$sample01$row))
  expect_equal(names(out_all$sample01$col), "foo")
  expect_equal(names(out_all$sample02$annot), "foo")
  expect_equal(out_all$sample01$col$foo, cgr1)
  expect_equal(out_all$sample02$annot$foo, agr2)
})

test_that("Get all applicable sample_ids in one margin", {
  out_mar1 <- spatialGraphs(sfe2, 1)
  expect_equal(names(out_mar1), c("sample01", "sample02"))
  expect_true(is.null(out_mar1$sample01))
  expect_true(is.null(out_mar1$sample02))

  out_mar2 <- spatialGraphs(sfe2, 2)
  expect_equal(names(out_mar2), c("sample01", "sample02"))
  expect_equal(names(out_mar2$sample01), "foo")
  expect_equal(out_mar2$sample02$bar, cgr2)

  out_mar3 <- spatialGraphs(sfe2, 3)
  expect_equal(names(out_mar3), c("sample01", "sample02"))
  expect_equal(names(out_mar3$sample01), "bar")
  expect_equal(out_mar3$sample02$foo, agr2)
})

test_that("Set all applicable margins in one sample_id", {
  out_sample02 <- spatialGraphs(sfe2, sample_id = "sample02")
  expect_equal(names(out_sample02), mar_names)
  expect_true(is.null(out_sample02$row))
  expect_equal(out_sample02$col$bar, cgr2)
  expect_equal(out_sample02$annot$foo, agr2)
})

test_that("Get all items in one margin and one sample_id", {
  out <- spatialGraphs(sfe2, 3, "sample01")
  expect_true(is.list(out))
  expect_equal(names(out), "bar")
  expect_equal(out$bar, agr1)
  # Make sure that it still works for one sample
  # I got into some troubles with this
  sfe <- sfe2[,colData(sfe2)$sample_id == "sample01"]
  out <- spatialGraphs(sfe2, 3, "sample01")
  expect_true(is.list(out))
  expect_equal(names(out), "bar")
  expect_equal(out$bar, agr1)
})

test_that("Get graphs of the same name from multiple samples", {
  colGraph(sfe2, "foo", sample_id = "sample02") <- cgr2
  out <- spatialGraphs(sfe2, 2, sample_id = "all", name = "foo")
  expect_true(is.list(out))
  expect_equal(length(out), 2L)
  expect_equal(names(out), c("sample01", "sample02"))
  expect_s3_class(out[[1]], "listw")
  expect_s3_class(out[[2]], "listw")
  # Still works with one sample
  sfe <- sfe2[,colData(sfe2)$sample_id == "sample01"]
  out <- spatialGraphs(sfe, 2, sample_id = "sample01", name = "foo")
  expect_true(is.list(out))
  expect_equal(length(out), 1L)
  expect_equal(names(out), "sample01")
  expect_s3_class(out[[1]], "listw")
})

test_that("Get one item in one margin and one sample_id", {
  out <- spatialGraph(sfe2, "foo", 2, "sample01")
  expect_equal(out, cgr1)
  # When the requested item doesn't exist
  expect_true(is.null(spatialGraph(sfe2, "bar", 1, "sample01")))
})

test_that("spatialGraphNames getter", {
  expect_equal(spatialGraphNames(sfe2, 2, "sample01"), "foo")
})

test_that("spatialGraphNames setter", {
  spatialGraphNames(sfe2, 2, "sample01") <- "meow"
  expect_equal(names(int_metadata(sfe2)$spatialGraphs["col", "sample01"][[1]]), "meow")
})

test_that("Set all graphs to NULL", {
  spatialGraphs(sfe2) <- NULL
  empty <- SpatialFeatureExperiment:::.initialize_spatialGraphs(sfe2)
  expect_equal(int_metadata(sfe2)$spatialGraphs, empty)
})

test_that("Set all graphs in one margin to NULL", {
  spatialGraphs(sfe2, 3) <- NULL
  expect_true(is.null(unlist(as.list(int_metadata(sfe2)$spatialGraphs[3,]))))
})

test_that("Set all graphs in one sample_id to NULL", {
  spatialGraphs(sfe2, sample_id = "sample01") <- NULL
  expect_true(is.null(unlist(as.list(int_metadata(sfe2)$spatialGraphs[,1]))))
})

test_that("Set all graphs in one margin of one sample_id to NULL", {
  spatialGraphs(sfe2, 3, "sample01") <- NULL
  expect_true(is.null(unlist(int_metadata(sfe2)$spatialGraphs[3,1])))
})

test_that("Remove one graph from one margin and one sample_id", {
  spatialGraph(sfe2, "foo", 3, "sample02") <- NULL
  expect_true(is.null(int_metadata(sfe2)$spatialGraphs[3,2][["foo"]]))
})
