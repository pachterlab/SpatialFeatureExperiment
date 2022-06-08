sfe_visium <- readRDS(system.file("testdata/sfe_visium.rds",
                                  package = "SpatialFeatureExperiment"))
g_visium <- readRDS(system.file("testdata/colgraph_visium.rds",
                                package = "SpatialFeatureExperiment"))
g_visium2 <- readRDS(system.file("testdata/colgraph_visium2.rds",
                                 package = "SpatialFeatureExperiment"))
colGraph(sfe_visium, "visium1", "sample01") <- g_visium
colGraph(sfe_visium, "visium2", "sample02") <- g_visium2
ag_samples <- readRDS(system.file("testdata/ag_samples.rds",
                                  package = "SpatialFeatureExperiment"))
annotGeometry(sfe_visium, "foo", sample_id = "all") <- ag_samples
sfe_visium1 <- sfe_visium[,colData(sfe_visium)$sample_id == "sample01"]
sfe_visium2 <- sfe_visium[,colData(sfe_visium)$sample_id == "sample02"]

sfe_concat <- cbind(sfe_visium1, sfe_visium2)

test_that("cbind properly deals with annotGeometries and spatialGraphs", {
  ag_orig <- annotGeometries(sfe_visium)
  ag_concat <- annotGeometries(sfe_concat)
  ag_concat <- ag_concat[names(ag_orig)]
  expect_equal(ag_orig, ag_concat)
  expect_equal(spatialGraphs(sfe_visium), spatialGraphs(sfe_concat))
})
