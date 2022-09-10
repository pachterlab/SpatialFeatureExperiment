# Unit test spdep and Visium graph wrappers
sfe2 <- readRDS(system.file("extdata/sfe_multi_sample.rds",
                            package = "SpatialFeatureExperiment"))
cgr1 <- readRDS(system.file("extdata/colgraph1.rds",
                            package = "SpatialFeatureExperiment"))

test_that("Get the correct graph and attr for reconstruction", {
  g <- findSpatialNeighbors(sfe2, "sample01", type = "spatialCoords",
                            MARGIN = 2, method = "tri2nb")
  expect_equal(g, cgr1, ignore_attr = TRUE)
  attrs_reconst <- attr(g, "method")
  expect_equal(names(attrs_reconst), c("FUN", "package", "args"))
  expect_equal(attrs_reconst$FUN, "findSpatialNeighbors")
  expect_equal(attrs_reconst$package, "SpatialFeatureExperiment")
  expect_true("row.names" %in% names(attrs_reconst$args))
  expect_equal(attrs_reconst$args$method, "tri2nb")
})

sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
                                  package = "SpatialFeatureExperiment"))
g_visium <- readRDS(system.file("extdata/colgraph_visium.rds",
                                package = "SpatialFeatureExperiment"))
test_that("Correct Visium graph", {
  g <- findVisiumGraph(sfe_visium, "sample01")
  expect_equal(g, g_visium, ignore_attr = TRUE)
  attrs_reconst <- attr(g, "method")
  expect_equal(attrs_reconst$FUN, "findVisiumGraph")
  expect_equal(attrs_reconst$args$style, "W")
})
