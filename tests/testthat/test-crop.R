library(sf)

sfe_visium <- readRDS(system.file("testdata/sfe_visium.rds",
                                  package = "SpatialFeatureExperiment"))
sfe_visium <- addVisiumSpotPoly(sfe_visium, 0.5)
bbox_use <- st_as_sfc(st_bbox(c(xmin = 1, xmax = 3, ymin = 1, ymax = 6), crs = NA))
bbox_use2 <- st_sf(geometry = bbox_use, sample_id = "sample01", crs = NA)

test_that("All spots in the cropped SFE objects indeed are covered by the bbox", {
  sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "all")
  expect_true(all(st_any_pred(spotPoly(sfe_cropped, "all"), bbox_use, pred = st_covered_by)))
  sfe_cropped2 <- crop(sfe_visium, y = bbox_use2, sample_id = "sample01")
  expect_true(all(st_any_pred(spotPoly(sfe_cropped2, "sample01"), bbox_use, pred = st_covered_by)))
  expect_false(all(st_any_pred(spotPoly(sfe_cropped2, "sample02"), bbox_use, pred = st_covered_by)))
  expect_equal(sum(st_any_intersects(spotPoly(sfe_cropped2, "sample02"), bbox_use)), 2)
})

annotGeometry(sfe_visium, "bbox", sample_id = "sample01") <- bbox_use2
test_that("annotPred", {
  out <- annotPred(sfe_visium, colGeometryName = "spotPoly",
                   annotGeometryName = "bbox", sample_id = "sample01")
  expect_equal(names(out), colnames(sfe_visium)[colData(sfe_visium)$sample_id == "sample01"])
  expect_true(all(out[c(1,2,5)]))
  expect_false(any(out[3:4]))
})

test_that("annotOp", {
  out <- annotOp(sfe_visium, colGeometryName = "spotPoly",
                 annotGeometryName = "bbox", sample_id = "sample01")
  expect_s3_class(out, "sf")
  expect_equal(rownames(out), colnames(sfe_visium)[colData(sfe_visium)$sample_id == "sample01"])
  p <- st_any_pred(out, bbox_use, st_covered_by)
  expect_true(all(p[c(1,2,5)]))
  expect_false(any(p[3:4]))
})

test_that("Find bbox of samples", {
  cg1 <- spotPoly(sfe_visium, "sample01")
  bbox1 <- bbox(sfe_visium, "sample01")
  expect_equal(st_bbox(st_union(cg1, bbox_use)), bbox1)
  bboxes <- bbox(sfe_visium, "all")
  expect_true(is.matrix(bboxes))
  expect_true(is.numeric(bboxes))
  expect_equal(colnames(bboxes), c("sample01", "sample02"))
  expect_equal(rownames(bboxes), c("xmin", "ymin", "xmax", "ymax"))
  expect_true(all(!is.na(bboxes)))
})

test_that("Remove empty space", {
  sfe_moved <- removeEmptySpace(sfe_visium, sample_id = "all")
  bboxes <- int_metadata(sfe_moved)$orig_bbox
  expect_true(is.matrix(bboxes))
  expect_true(is.numeric(bboxes))
  expect_equal(colnames(bboxes), c("sample01", "sample02"))
  expect_equal(rownames(bboxes), c("xmin", "ymin", "xmax", "ymax"))
  expect_true(all(!is.na(bboxes)))
  new_bboxes <- bbox(sfe_moved, "all")
  expect_true(all(abs(new_bboxes[c("xmin", "ymin"), c("sample01", "sample02")]) < .Machine$double.eps))
})
