library(sf)

sfe_visium <- readRDS(system.file("testdata/sfe_visium.rds",
                                  package = "SpatialFeatureExperiment"))
sfe_visium <- addVisiumSpotPoly(sfe_visium, 0.5)

test_that("All spots in the cropped SFE objects indeed intersects the bbox", {
  bbox_use <- st_as_sfc(st_bbox(c(xmin = 1, xmax = 3, ymin = 1, ymax = 6), crs = NA))
  sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "all")
  expect_true(all(st_any_pred(spotPoly(sfe_cropped, "all"), bbox_use, pred = st_covered_by)))
  bbox_use2 <- st_sf(geometry = bbox_use, sample_id = "sample01", crs = NA)
  sfe_cropped2 <- crop(sfe_visium, y = bbox_use2, sample_id = "sample01")
  expect_true(all(st_any_pred(spotPoly(sfe_cropped2, "sample01"), bbox_use, pred = st_covered_by)))
  expect_false(all(st_any_pred(spotPoly(sfe_cropped2, "sample02"), bbox_use, pred = st_covered_by)))
  expect_equal(sum(st_any_intersects(spotPoly(sfe_cropped2, "sample02"), bbox_use)), 2)
})

test_that("Find bbox of samples", {
  cg1 <- spotPoly(sfe_visium, "sample01")
  bbox1 <- bbox(sfe_visium, "sample01")
  expect_equal(st_bbox(cg1), bbox1)
  bboxes <- bbox(sfe_visium, "all")
  expect_s3_class(bboxes, "matrix")
  expect_equal(colnames(bboxes), c("sample01", "sample02"))
  expect_equal(rownames(bboxes), c("xmin", "ymin", "xmax", "ymax"))
})

test_that("Remove empty space", {
  sfe_moved <- removeEmptySpace(sfe_visium, sample_id = "all")

})
