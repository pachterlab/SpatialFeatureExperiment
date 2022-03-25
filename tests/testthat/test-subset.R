# Unit test the subsetting method
# 4. Warning message and dropping graphs when reconstruction info is unavailable
# 5. Warning message and dropping graphs when package required for reconstruction is not installed
sfe2 <- readRDS(system.file("testdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))
agr1 <- readRDS(system.file("testdata/annotgraph1.rds",
                            package = "SpatialFeatureExperiment"))
ag <- readRDS(system.file("testdata/ag.rds",
                          package = "SpatialFeatureExperiment"))
# Should have passed unit test for graph_wrapper
spatialGraph(sfe2, type = "foo", MARGIN = 2, sample_id = "sample01") <-
  findSpatialNeighbors(sfe2, "sample01", type = "spatialCoords", MARGIN = 2,
                       method = "tri2nb")
annotGraph(sfe2, "bar", "sample01") <- agr1
annotGeometry(sfe2, "baz", "sample01") <- ag

test_that("After removing one sample_id, it's also removed in annotGeometries", {
  sfe2 <- sfe2[,4:5]
  expect_true(is.null(int_metadata(sfe2)$annotGeometries$baz))
})

test_that("row and col graphs are dropped if reconstruct_graphs = FALSE", {
  expect_warning(sfe2 <- sfe2[, 2:5, reconstruct_graphs = FALSE], "Dropping all")
  # Don't have rowGraphs to begin with
  expect_true(is.null(unlist(as.list(colGraphs(sfe2)))))
})

test_that("Correctly reconstruct the graphs when they need to be reconstructed", {

})
