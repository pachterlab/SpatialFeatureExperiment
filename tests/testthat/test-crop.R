library(sf)

sfe_visium <- readRDS(system.file("testdata/sfe_visium.rds",
                                  package = "SpatialFeatureExperiment"))
sfe_visium <- addVisiumSpotPoly(sfe_visium, 0.5)
bbox_use <- st_as_sfc(st_bbox(c(xmin = 1, xmax = 3, ymin = 1, ymax = 6), crs = NA))
sfe_cropped <- crop(sfe_visium, bbox_use)
test_that("All spots in the cropped SFE objects indeed intersects the bbox", {
  expect_true(all(st_any_intersects(spotPoly(sfe_cropped, "all"), bbox_use)))
})
