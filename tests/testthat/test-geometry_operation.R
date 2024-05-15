library(sf)
library(SFEData)
# From sf's examples
pts = st_sfc(st_point(c(.5,.5)), st_point(c(1.5, 1.5)), st_point(c(2.5, 2.5)))
pol = st_polygon(list(rbind(c(0,0), c(2,0), c(2,2), c(0,2), c(0,0))))

test_that("trivial, st_any_intersects", {
    o <- st_any_intersects(pts, pol)
    expect_equal(o, c(TRUE, TRUE, FALSE))
})

test_that("trivial, st_n_intersects", {
    n1 <- st_n_intersects(pts, pol)
    expect_equal(n1, c(1,1,0))
    n2 <- st_n_intersects(pol, pts)
    expect_equal(n2, 2)
})

# Crop============

sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
    package = "SpatialFeatureExperiment"
))
sfe_visium <- addVisiumSpotPoly(sfe_visium, 0.5)
bbox_use <- st_as_sfc(st_bbox(c(xmin = 1, xmax = 3, ymin = 1, ymax = 6),
    crs = NA
))
bbox_use2 <- st_as_sfc(st_bbox(c(xmin = 1, xmax = 3, ymin = 5, ymax = 9)), crs = NA)
bbox_sf <- st_sf(geometry = bbox_use, sample_id = "sample01", crs = NA)

cg <- colGeometry(sfe_visium, sample_id = "sample01")
bbox_cg <- st_bbox(st_centroid(cg))
ag <- st_bbox(cg) |> st_as_sfc() |> st_buffer(dist = 0.3)
annotGeometry(sfe_visium, "box", sample_id = "sample01") <- st_sf(geometry = ag, sample_id = "sample01",
                                                                  sf_column_name = "geometry",
                                                                  crs = NA)
rg_bbox <- st_bbox(ag)
set.seed(29)
rg1 <- cbind(runif(20, rg_bbox["xmin"], rg_bbox["xmax"]),
             runif(20, rg_bbox["ymin"], rg_bbox["ymax"]))
rg2 <- matrix(rnorm(30, mean = 3, sd = 3), ncol = 2)
rg_use <- st_sf(geometry = st_sfc(st_multipoint(rg1), st_multipoint(rg2)),
                crs = NA)
rowGeometry(sfe_visium, "points", sample_id = "sample01", withDimnames = FALSE) <- rg_use

test_that("All spots in the cropped SFE objects indeed are covered by the bbox", {
    sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "all")
    cg <- spotPoly(sfe_cropped, "all")
    expect_true(all(st_any_pred(cg, bbox_use, pred = st_covered_by)))
    expect_true(st_geometry_type(cg, by_geometry = FALSE) == "POLYGON")
    rg <- rowGeometry(sfe_cropped)
    expect_true(all(st_any_pred(rg, bbox_use, pred = st_covered_by)))
})

test_that("Only crop one sample out of two, with sf", {
    expect_error(crop(sfe_visium, y = bbox_sf, sample_id = "sample02"),
                 "No bounding boxes for samples specified.")
    sfe_cropped2 <- crop(sfe_visium, y = bbox_sf, sample_id = "sample01")
    expect_true(all(st_any_pred(spotPoly(sfe_cropped2, "sample01"), bbox_use,
                                pred = st_covered_by
    )))
    expect_false(all(st_any_pred(spotPoly(sfe_cropped2, "sample02"), bbox_use,
                                 pred = st_covered_by
    )))
    expect_equal(sum(st_any_intersects(
        spotPoly(sfe_cropped2, "sample02"),
        bbox_use
    )), 2)
})

test_that("Using a bounding box to crop SFE objects, deprecated way", {
    expect_warning(sfe_cropped <- crop(sfe_visium, sample_id = "sample01",
                                       xmin = 1, xmax = 3, ymin = 1, ymax = 6),
                   "deprecated")
    cg <- spotPoly(sfe_cropped, "sample01")
    expect_true(all(st_any_pred(cg, bbox_use, pred = st_covered_by)))
})

test_that("Using a bounding box to crop SFE objects, current way, expected errors", {
    expect_error(crop(sfe_visium, y = "foobar", sample_id = "sample01"),
                 "bbox must be a numeric vector or matrix.")
    expect_error(crop(sfe_visium, y = c(meow = 1, purr = 2), sample_id = "sample01"),
                 "must be a vector of length 4")
    m <- matrix(1:8, ncol = 2)
    expect_error(crop(sfe_visium, y = m, sample_id = "sample01"),
                 "must have rownames xmin, xmax")
    rownames(m) <- c("xmin", "ymin", "xmax", "ymax")
    expect_error(crop(sfe_visium, y = m, sample_id = "sample01"),
                 "must have colnames")
})

test_that("Using a bounding box to crop SFE objects, current way, one sample", {
    # Sample ID for bbox doesn't match sample specified in crop function
    m <- matrix(c(1, 3, 1, 6), ncol = 1,
                dimnames = list(c("xmin", "xmax", "ymin", "ymax"),
                                "sample01"))
    expect_error(crop(sfe_visium, y = m, sample_id = "sample02"),
                 "No bounding boxes for samples specified.")
    sfe_cropped <- crop(sfe_visium, y = m)
    expect_true(all(st_any_pred(spotPoly(sfe_cropped, "sample01"), bbox_use,
                                pred = st_covered_by
    )))
    expect_false(all(st_any_pred(spotPoly(sfe_cropped, "sample02"), bbox_use,
                                 pred = st_covered_by
    )))
    expect_equal(sum(st_any_intersects(
        spotPoly(sfe_cropped, "sample02"),
        bbox_use
    )), 2)
})

# I don't expect it to be common to use the same bbox to crop all samples,
# because the ROIs in different samples most likely have different bboxes
test_that("Using a bounding box to crop SFE objects, current way, all samples", {
    m <- matrix(c(1, 3, 1, 6, 1, 3, 5, 9), ncol = 2,
                dimnames = list(c("xmin", "xmax", "ymin", "ymax"),
                                c("sample01", "sample02")))
    sfe_cropped <- crop(sfe_visium, y = m)
    expect_true(all(st_any_pred(spotPoly(sfe_cropped, "sample01"), bbox_use,
                                pred = st_covered_by
    )))
    expect_true(all(st_any_pred(spotPoly(sfe_cropped, "sample02"), bbox_use2,
                                pred = st_covered_by
    )))
})

test_that("When a geometry is broken into multiple pieces", {
    notch <- st_as_sfc(st_bbox(c(xmin = 0, xmax = 1.5, ymin = 1.7, ymax = 1.9)))
    bbox_use3 <- st_difference(bbox_use, notch)
    sfe_cropped3 <- crop(sfe_visium, bbox_use3, sample_id = "sample01")
    cg <- spotPoly(sfe_cropped3, "all")
    expect_true(st_geometry_type(cg, by_geometry = FALSE) == "MULTIPOLYGON")
})

test_that("When a sample is removed by cropping", {
    m <- matrix(c(1, 3, 1, 6, 6, 9, 10, 13), ncol = 2,
                dimnames = list(c("xmin", "xmax", "ymin", "ymax"),
                                c("sample01", "sample02")))
    expect_warning(sfe_cropped <- crop(sfe_visium, m), "were removed")
    expect_equal(sampleIDs(sfe_cropped), "sample01")
})

test_that("Keep whole colGeometry items", {
    # In case you don't want small slivers or broken into multiple pieces
    # rowGeometries, annotGeometries, and images should be cropped by the actual
    # bbox of the remaining items of the colGeometry
    sfe_cropped <- crop(sfe_visium, bbox_cg, sample_id = "sample01",
                        keep_whole = "col")
    cg2 <- colGeometry(sfe_cropped, sample_id = "sample01")
    # colGeometry not cropped
    expect_equal(cg, cg2)
    # annotGeometry cropped by bbox of colGeometry
    cg_bbox <- st_bbox(cg) |> st_as_sfc()
    ag2 <- annotGeometry(sfe_cropped, sample_id = "sample01")
    expect_true(st_equals(ag2, cg_bbox, sparse = FALSE))
    # rowGeometry
    rg_check <- rowGeometry(sfe_cropped, "points", sample_id = "sample01")
    rg_bbox2 <- st_bbox(rg_check) |> st_as_sfc()
    expect_true(st_covered_by(rg_bbox2, cg_bbox, sparse = FALSE))
})

test_that("Keep whole annotGeometry items", {
    cg <- colGeometry(sfe_visium, sample_id = "sample01")
    bbox_use <- st_bbox(st_centroid(cg))
    set.seed(29)
    ag_use <- cbind(runif(20, bbox_use["xmin"], bbox_use["xmax"]),
                    runif(20, bbox_use["ymin"], bbox_use["ymax"])) |>
        as.data.frame()
    names(ag_use) <- c("x", "y")
    ag_use <- st_as_sf(ag_use, coords = c("x", "y"), crs = NA)
    ag_use <- st_buffer(ag_use, 0.5)
    ag_use$sample_id <- "sample01"
    annotGeometry(sfe_visium, "circles", "sample01") <- ag_use
    sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "sample01",
                        keep_whole = "annot")
    expect_true(st_covered_by(spotPoly(sfe_cropped, "sample01"), st_as_sfc(bbox_use),
                              sparse = FALSE) |> all())
    expect_equal(ag_use, annotGeometry(sfe_cropped, "circles", "sample01")[,names(ag_use)])
})

test_that("Use st_difference for cropping, not cover", {
    bbox_use <- c(xmin = 2.5, xmax = 3.5, ymin = 1.75, ymax = 2.5)
    cg <- colGeometry(sfe_visium, sample_id = "sample01")
    bbox1 <- st_bbox(st_centroid(cg))
    set.seed(29)
    ag_use <- cbind(runif(20, bbox1["xmin"], bbox1["xmax"]),
                    runif(20, bbox1["ymin"], bbox1["ymax"])) |>
        as.data.frame()
    names(ag_use) <- c("x", "y")
    ag_use <- st_as_sf(ag_use, coords = c("x", "y"), crs = NA)
    ag_use <- st_buffer(ag_use, 0.5)
    ag_use$sample_id <- "sample01"
    annotGeometry(sfe_visium, "circles", "sample01") <- ag_use
    sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "sample01",
                        op = st_difference)
    cg2 <- spotPoly(sfe_cropped, "sample01")
    ag2 <- annotGeometry(sfe_cropped, "circles", "sample01")
    expect_false(any(st_overlaps(cg2, st_as_sfc(st_bbox(bbox_use)), sparse = FALSE)))
    expect_false(any(st_covered_by(cg2, st_as_sfc(st_bbox(bbox_use)), sparse = FALSE)))
    expect_false(any(st_overlaps(ag2, st_as_sfc(st_bbox(bbox_use)), sparse = FALSE)))
    expect_false(any(st_covered_by(ag2, st_as_sfc(st_bbox(bbox_use)), sparse = FALSE)))
})

test_that("Use st_difference for cropping, cover", {
    bbox_use <- c(xmin = 2.5, xmax = 3.5, ymin = 1.75, ymax = 2.5)
    cg <- colGeometry(sfe_visium, sample_id = "sample01")
    bbox1 <- st_bbox(st_centroid(cg)) |> st_as_sfc()
    set.seed(29)
    ag_use <- st_sample(bbox1, 20)
    ag_use <- st_sf(geometry = ag_use, coords = c("x", "y"), crs = NA)
    ag_use <- st_buffer(ag_use, 0.5)
    ag_use$sample_id <- "sample01"
    annotGeometry(sfe_visium, "circles", "sample01") <- ag_use
    sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "sample01",
                        op = st_difference, cover = TRUE,
                        keep_whole = c("col", "annot"))
    cg2 <- spotPoly(sfe_cropped, "sample01")
    ag2 <- annotGeometry(sfe_cropped, "circles", "sample01")
    bbox_sf <- st_as_sfc(st_bbox(bbox_use))
    expect_true(all(st_disjoint(cg2, bbox_sf, sparse = FALSE)))
    expect_true(all(st_disjoint(ag2, bbox_sf, sparse = FALSE)))
})

test_that("Only cells/spots covered by y if keep whole", {
    cg <- colGeometry(sfe_visium, sample_id = "sample01")
    bbox1 <- st_bbox(st_centroid(cg))
    bbox2 <- st_bbox(cg)
    bbox_use <- c(bbox2[c("xmin", "xmax", "ymin")], bbox1["ymax"])
    set.seed(29)
    ag_use <- cbind(runif(20, bbox_use["xmin"], bbox_use["xmax"]),
                    runif(20, bbox_use["ymin"], bbox_use["ymax"])) |>
        as.data.frame()
    names(ag_use) <- c("x", "y")
    ag_use <- st_as_sf(ag_use, coords = c("x", "y"), crs = NA)
    ag_use <- st_buffer(ag_use, 0.5)
    ag_use$sample_id <- "sample01"
    annotGeometry(sfe_visium, "circles", "sample01") <- ag_use
    sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "sample01",
                        keep_whole = c("col", "annot"), cover = TRUE)
    cg2 <- spotPoly(sfe_cropped, "sample01")
    expect_true(all(st_covered_by(cg2, st_as_sfc(st_bbox(bbox_use)), sparse = FALSE)))
    expect_true(all(lengths(st_equals(cg2, cg))))
    ag2 <- annotGeometry(sfe_cropped, "circles", "sample01")
    expect_true(all(st_covered_by(ag2, st_as_sfc(st_bbox(bbox_use)), sparse = FALSE)))
    expect_true(all(lengths(st_equals(ag2, ag_use))))
})

test_that("Error when other spatial operations are specified", {
    bbox_use <- c(xmin = 2.5, xmax = 3.5, ymin = 1.75, ymax = 2.5)
    expect_error(crop(sfe_visium, bbox_use, sample_id = "sample01",
                      op = st_sym_difference),
                 "op must be either st_intersection or st_difference")
})

test_that("Crop 3D geometry", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))
    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     z_option = "3d")
    bbox1 <- c(xmin = 171500, ymin = 11500, xmax = 172000, ymax = 12000)
    sfe_cropped <- crop(sfe, bbox1)
    bbox_new <- bbox(sfe_cropped)
    expect_true(st_covered_by(st_as_sfc(st_bbox(bbox_new)), st_as_sfc(st_bbox(bbox1)),
                              sparse = FALSE))
    rg <- txSpots(sfe_cropped)
    expect_true(st_covered_by(st_as_sfc(st_bbox(rg)), st_as_sfc(st_bbox(bbox_new)),
                              sparse = FALSE))
    expect_equal(unclass(st_z_range(rg)), c(zmin = 0, zmax = 1),
                 ignore_attr = "crs")
    unlink(dir_use, recursive = TRUE)
})

# annotPred and annotOp===========
annotGeometry(sfe_visium, "bbox", sample_id = "sample01") <- bbox_sf
test_that("annotPred", {
    out <- annotPred(sfe_visium,
        colGeometryName = "spotPoly",
        annotGeometryName = "bbox", sample_id = "sample01"
    )
    expect_equal(
        names(out),
        colnames(sfe_visium)[colData(sfe_visium)$sample_id == "sample01"]
    )
    expect_true(all(out[c(1, 2, 5)]))
    expect_false(any(out[3:4]))
})

test_that("annotOp", {
    out <- annotOp(sfe_visium,
        colGeometryName = "spotPoly",
        annotGeometryName = "bbox", sample_id = "sample01"
    )
    expect_s3_class(out, "sf")
    expect_equal(
        rownames(out),
        colnames(sfe_visium)[colData(sfe_visium)$sample_id == "sample01"]
    )
    p <- st_any_pred(out, bbox_use, st_covered_by)
    expect_true(all(p[c(1, 2, 5)]))
    expect_false(any(p[3:4]))
})

# bbox=============
test_that("Find bbox of samples", {
    cg1 <- spotPoly(sfe_visium, "sample01")
    bbox1 <- bbox(sfe_visium, "sample01")
    expect_equal(st_bbox(rowGeometry(sfe_visium)), bbox1,
        ignore_attr = c("class", "crs")
    )
    bboxes <- bbox(sfe_visium, "all")
    expect_true(is.matrix(bboxes))
    expect_true(is.numeric(bboxes))
    expect_equal(colnames(bboxes), c("sample01", "sample02"))
    expect_equal(rownames(bboxes), c("xmin", "ymin", "xmax", "ymax"))
    expect_true(all(!is.na(bboxes)))
})

test_that("bbox when 0 rows or columns", {
    sfe0 <- sfe_visium[,logical(0)]
    bbox0 <- bbox(sfe0)
    expect_equal(names(bbox0), c("xmin", "ymin", "xmax", "ymax"))
    expect_true(all(is.na(bbox0)))
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
    expect_true(all(abs(new_bboxes[
        c("xmin", "ymin"),
        c("sample01", "sample02")
    ]) < .Machine$double.eps))
    # Make sure that spatialCoords are correctly moved
    coord_diffs <- spatialCoords(sfe_visium) - spatialCoords(sfe_moved)
    diff1 <- unname(coord_diffs[sfe_visium$sample_id == "sample01",][1,])
    diff2 <- unname(coord_diffs[sfe_visium$sample_id == "sample02",][1,])
    expect_equal(diff1, unname(bboxes[c("xmin", "ymin"), "sample01"]))
    expect_equal(diff2, unname(bboxes[c("xmin", "ymin"), "sample02"]))
})

# removeEmptySpace===============
library(SFEData)
sfe1 <- McKellarMuscleData("small")
sfe2 <- McKellarMuscleData("small2")
sfe <- cbind(sfe1, sfe2)

# Toy example before removing empty space
bboxes <- bbox(sfe, sample_id = "all")
set.seed(29)
cg_toy1 <- data.frame(
    x = runif(ncol(sfe1), bboxes["xmin", "Vis5A"], bboxes["xmax", "Vis5A"]),
    y = runif(ncol(sfe1), bboxes["ymin", "Vis5A"], bboxes["ymax", "Vis5A"]),
    sample_id = "Vis5A")
cg_toy1 <- df2sf(cg_toy1)
cg_toy2 <- data.frame(
    x = runif(ncol(sfe2), bboxes["xmin", "sample02"], bboxes["xmax", "sample02"]),
    y = runif(ncol(sfe2), bboxes["ymin", "sample02"], bboxes["ymax", "sample02"]),
    sample_id = "sample02")
cg_toy2 <- df2sf(cg_toy2)
cg_toy <- rbind(cg_toy1, cg_toy2)
sfe <- sfe[rowSums(counts(sfe)) > 0,]
sfe_shifted <- removeEmptySpace(sfe)
bbox_new <- st_as_sfc(st_bbox(bbox(sfe_shifted, sample_id = "Vis5A")))
bbox_new2 <- st_as_sfc(st_bbox(bbox(sfe_shifted, sample_id = "sample02")))
bbox_old1 <- st_as_sfc(st_bbox(bbox(sfe, sample_id = "Vis5A")))
bbox_old2 <- st_as_sfc(st_bbox(bbox(sfe, sample_id = "sample02")))

test_that("colGeometry setter after removing empty space", {
    # Two samples
    colGeometry(sfe_shifted, "toy", sample_id = "all",
                withDimnames = FALSE) <- cg_toy
    cg <- colGeometry(sfe_shifted, "toy", sample_id = "all")
    bb <- st_as_sfc(st_bbox(cg))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # One of two samples
    colGeometry(sfe_shifted, "toy2", sample_id = "sample02",
                withDimnames = FALSE) <- cg_toy2
    cg <- colGeometry(sfe_shifted, "toy2", sample_id = "sample02")
    bb <- st_as_sfc(st_bbox(cg))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # Use colGeometries setter
    colGeometries(sfe_shifted, withDimnames = FALSE) <- list(foo = cg_toy)
    cg <- colGeometry(sfe_shifted, "foo", sample_id = "all")
    bb <- st_as_sfc(st_bbox(cg))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # Not moved with translate = FALSE
    colGeometry(sfe_shifted, "bar", sample_id = "all", translate = FALSE,
                withDimnames = FALSE) <- cg_toy
    cg1 <- colGeometry(sfe_shifted, "bar", sample_id = "Vis5A")
    cg2 <- colGeometry(sfe_shifted, "bar", sample_id = "sample02")
    bb1 <- st_as_sfc(st_bbox(cg1))
    bb2 <- st_as_sfc(st_bbox(cg2))
    expect_true(st_covered_by(bb1, bbox_old1, sparse = FALSE))
    expect_true(st_covered_by(bb2, bbox_old2, sparse = FALSE))
})

test_that("annotGeometry setter after removing empty space", {
    # One sample, that it shifted
    sfe_shifted1 <- removeEmptySpace(sfe1)
    ag <- annotGeometry(sfe_shifted1, "myofiber_simplified")
    bb <- st_as_sfc(st_bbox(ag))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # One sample setter
    ag_old <- annotGeometry(sfe1, "myofiber_simplified")
    annotGeometry(sfe_shifted1, "foo") <- ag_old
    ag <- annotGeometry(sfe_shifted1, "foo")
    bb <- st_as_sfc(st_bbox(ag))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # Two samples
    ag_old2 <- annotGeometry(sfe, "myofiber_simplified", sample_id = "all")
    annotGeometry(sfe_shifted, "foo", sample_id = "all") <- ag_old2
    ag <- annotGeometry(sfe_shifted, "foo", sample_id = "all")
    bb <- st_as_sfc(st_bbox(ag))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # One of two samples to modify existing annotGeometry
    annotGeometry(sfe_shifted, "bar", sample_id = "Vis5A") <- ag_old
    ag <- annotGeometry(sfe_shifted, "foo", sample_id = "Vis5A")
    bb <- st_as_sfc(st_bbox(ag))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # Use annotGeometries setter
    annotGeometries(sfe_shifted) <- list(foobar = ag_old2)
    ag <- annotGeometry(sfe_shifted, "foobar", sample_id = "all")
    bb <- st_as_sfc(st_bbox(ag))
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))

    # Not moved when translate = FALSE
    annotGeometry(sfe_shifted, "baz", sample_id = "all", translate = FALSE) <-
        ag_old2
    bb1 <- st_as_sfc(st_bbox(annotGeometry(sfe_shifted, "baz", sample_id = "Vis5A")))
    bb2 <- st_as_sfc(st_bbox(annotGeometry(sfe_shifted, "baz", sample_id = "sample02")))
    expect_true(st_covered_by(bb1, bbox_old1, sparse = FALSE))
    expect_true(st_covered_by(bb2, bbox_old2, sparse = FALSE))
})

set.seed(29)
rg1 <- st_sf(geometry = st_sample(bbox_old1, nrow(sfe)))
rg2 <- st_sf(geometry = st_sample(bbox_old2, nrow(sfe)))

rowGeometry(sfe, "foo", "Vis5A", withDimnames = FALSE) <- rg1
rowGeometry(sfe, "foo", "sample02", withDimnames = FALSE) <- rg2

test_that("rowGeometry removing empty space", {
    sfe_shifted <- removeEmptySpace(sfe)
    rg1_sh <- rowGeometry(sfe_shifted, "foo", "Vis5A")
    rg2_sh <- rowGeometry(sfe_shifted, "foo", "sample02")
    bb1 <- st_bbox(rg1_sh) |> st_as_sfc()
    bb2 <- st_bbox(rg2_sh) |> st_as_sfc()
    expect_true(st_covered_by(bb1, bbox_new, sparse = FALSE))
    expect_true(st_covered_by(bb2, bbox_new2, sparse = FALSE))
})

test_that("rowGeometry setter after removing empty space", {
    rowGeometry(sfe_shifted, "bar", "Vis5A", withDimnames = FALSE) <- rg1
    rg_shifted <- rowGeometry(sfe_shifted, "bar", "Vis5A")
    bb <- st_bbox(rg_shifted) |> st_as_sfc()
    expect_true(st_covered_by(bb, bbox_new, sparse = FALSE))
})

test_that("Don't translate value if it's already translated", {
    bb <- st_bbox(cg_toy)
    cg_translated <- cg_toy
    cg_translated$geometry <- cg_translated$geometry - bb[c("xmin", "ymin")]
    bb2 <- st_as_sfc(st_bbox(cg_translated))
    colGeometry(sfe_shifted, "cg_translated", sample_id = "all",
                translate = TRUE, withDimnames = FALSE) <- cg_translated
    cg_check <- colGeometry(sfe_shifted, "cg_translated", "all")
    expect_true(all(st_covered_by(cg_check, bb2, sparse = FALSE)))
})

# annotSummary========
test_that("annotSummary", {
    out <- annotSummary(sfe1, "spotPoly", "myofiber_simplified", "area")
    expect_s3_class(out, "data.frame")
    expect_equal(rownames(out), colnames(sfe1))
    expect_equal(names(out), "area")
    expect_true(is.numeric(out$area))
})

# Operations when there're images=================
# Need uncropped image
if (!dir.exists("ob")) dir.create(file.path("ob", "outs"), recursive = TRUE)
mat_fn <- file.path("ob", "outs", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
                  destfile = file.path("ob", "outs", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("ob", "outs", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz",
                  destfile = file.path("ob", "outs", "spatial.tar.gz"))
    untar(file.path("ob", "outs", "spatial.tar.gz"), exdir = file.path("ob", "outs"))
}

library(SpatialExperiment)
library(terra)
library(SingleCellExperiment)
library(S4Vectors)
sfe <- read10xVisiumSFE("ob")

test_that("bbox when images are included", {
    bbox_tot <- bbox(sfe, include_image = TRUE) |> st_bbox() |> st_as_sfc()
    bbox_img <- ext(getImg(sfe, image_id = "hires")) |>
        st_bbox() |> st_as_sfc()
    expect_true(st_equals(bbox_tot, bbox_img, sparse = FALSE))
})

sfe <- sfe[rowSums(counts(sfe)) > 0,]
bbox_use <- bbox(sfe) |> st_bbox() |> st_as_sfc()
set.seed(29)
ag <- st_sample(bbox_use, 20) |> st_buffer(dist = 100)
ag <- st_sf(geometry = ag, sample_id = "sample01")
rg <- st_sample(bbox_use, nrow(sfe))
rg <- st_sf(geometry = rg)
rownames(rg) <- rownames(sfe)
annotGeometry(sfe, "foo") <- ag
rowGeometry(sfe, "bar") <- rg

test_that("When no cells/spots left after cropping", {
    bbox0 <- c(xmin = 3000, xmax = 4000, ymin = 2000, ymax = 3000)
    expect_warning(sfe0 <- SpatialFeatureExperiment::crop(sfe, bbox0),
                   "No cells/spots left after cropping")
    expect_equal(ncol(sfe0), 0)
    expect_equal(nrow(sfe0), nrow(sfe))
    expect_equal(nrow(spotPoly(sfe0)), 0)
    expect_true(isEmpty(rowGeometries(sfe0)))
    expect_equal(nrow(imgData(sfe0)), 0)
    expect_equal(annotGeometryNames(sfe0), annotGeometryNames(sfe))
    expect_equal(nrow(annotGeometry(sfe0, "foo")), 0)
})

test_that("Image is shifted after removing empty space", {
    sfe2 <- removeEmptySpace(sfe)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
    expect_true(st_area(bbox_geom) / st_area(bbox_img) > 0.97)
})

test_that("Image is cropped after cropping SFE object", {
    bc <- bbox(sfe)
    bbox_use <- c(xmin = bc["xmin"], xmax = bc["xmin"] + 2000,
                  ymin = bc["ymin"], ymax = bc["ymin"] + 2000) |>
        setNames(c("xmin", "xmax", "ymin", "ymax")) |>
        st_bbox() |> st_as_sfc()
    sfe2 <- SpatialFeatureExperiment::crop(sfe, bbox_use)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2, image_id = "hires"))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean, use = "complete.obs")) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
    expect_true(st_area(bbox_geom) / st_area(bbox_img) > 0.99)
})

# Affine transformations of the SFE object==============
test_that("Transpose SFE object with image", {
    sfe2 <- transpose(sfe)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_cg <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(rowGeometry(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Transpose SFE object with image, after cropping image", {
    sfe <- sfe[,sfe$in_tissue]
    sfe2 <- transpose(sfe)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2, image_id = "hires"))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    bbox_cg_orig <- st_bbox(spotPoly(sfe))
    bbox_cg_orig_sf <- bbox_cg_orig |> st_as_sfc()
    bbox_img_orig <- ext(getImg(sfe, image_id = "hires"))
    bbox_img_orig_sf <- bbox_img_orig |> st_bbox() |> st_as_sfc()
    bbox_cg <- st_bbox(spotPoly(sfe2))
    bbox_cg_sf <- bbox_cg |> st_as_sfc()
    bbox_img <- as.vector(ext(img))
    bbox_img_sf <- bbox_img |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_cg_sf, bbox_img_sf, sparse = FALSE))
    expect_equal(bbox_img_orig[["ymax"]] - bbox_cg_orig[["ymax"]],
                 bbox_cg[["xmin"]] - bbox_img[["xmin"]])

    bbox_rg <- st_bbox(rowGeometry(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg_sf, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with image, vertical", {
    sfe2 <- mirror(sfe, direction = "vertical")
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_cg <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(rowGeometry(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with image after cropping", {
    sfe <- sfe[,sfe$in_tissue]
    sfe2 <- mirror(sfe, direction = "vertical")
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2, image_id = "hires"))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    bbox_cg_orig <- st_bbox(spotPoly(sfe))
    bbox_cg_orig_sf <- bbox_cg_orig |> st_as_sfc()
    bbox_img_orig <- ext(getImg(sfe, image_id = "hires"))
    bbox_img_orig_sf <- bbox_img_orig |> st_bbox() |> st_as_sfc()
    bbox_cg <- st_bbox(spotPoly(sfe2))
    bbox_cg_sf <- bbox_cg |> st_as_sfc()
    bbox_img <- as.vector(ext(img))
    bbox_img_sf <- bbox_img |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_cg_sf, bbox_img_sf, sparse = FALSE))
    expect_equal(bbox_img_orig[["ymax"]] - bbox_cg_orig[["ymax"]],
                 st_bbox(bbox_cg_sf)[["ymin"]] - bbox_img[["ymin"]])

    bbox_rg <- st_bbox(rowGeometry(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg_sf, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with image, horizontal", {
    sfe2 <- mirror(sfe, direction = "horizontal")
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_cg <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(rowGeometry(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Rotate SFE object with image", {
    sfe2 <- SpatialFeatureExperiment::rotate(sfe, degrees = 45)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2) |> toSpatRasterImage(save_geotiff = FALSE) |> imgRaster()

    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    # Also see if spotPoly is aligned
    v2 <- terra::extract(terra::mean(img), st_centroid(spotPoly(sfe2)))
    expect_true(abs(cor(nCounts, v2$mean)) > 0.4)

    int1 <- st_intersects(rg, spotPoly(sfe))
    int2 <- st_intersects(rowGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)

    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Rotate SFE object with image after cropping", {
    sfe <- sfe[,sfe$in_tissue]
    sfe2 <- SpatialFeatureExperiment::rotate(sfe, degrees = 45)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2) |> toSpatRasterImage(save_geotiff = FALSE) |> imgRaster()

    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    # Also see if spotPoly is aligned
    v2 <- terra::extract(terra::mean(img), st_centroid(spotPoly(sfe2)))
    expect_true(abs(cor(nCounts, v2$mean)) > 0.4)

    int1 <- st_intersects(rg, spotPoly(sfe))
    int2 <- st_intersects(rowGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)

    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Scale SFE object with image", {
    sfe2 <- SpatialFeatureExperiment::scale(sfe, factor = 1.5)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2) |> imgRaster()

    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    # Also see if spotPoly is aligned
    v2 <- terra::extract(terra::mean(img), st_centroid(spotPoly(sfe2)))
    expect_true(abs(cor(nCounts, v2$mean)) > 0.4)

    int1 <- st_intersects(rg, spotPoly(sfe))
    int2 <- st_intersects(rowGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)

    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Scale SFE object with image after cropping", {
    sfe <- sfe[,sfe$in_tissue]
    sfe2 <- SpatialFeatureExperiment::scale(sfe, factor = 1.5)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2) |> imgRaster()

    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    # Also see if spotPoly is aligned
    v2 <- terra::extract(terra::mean(img), st_centroid(spotPoly(sfe2)))
    expect_true(abs(cor(nCounts, v2$mean)) > 0.4)

    int1 <- st_intersects(rg, spotPoly(sfe))
    int2 <- st_intersects(rowGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)

    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("General affine transformation of SFE object with image", {
    M <- matrix(c(0.6, -0.2, 0.2, 0.3), nrow = 2)
    v <- c(0, 300)
    sfe2 <- SpatialFeatureExperiment::affine(sfe, M = M, v = v)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2) |> toSpatRasterImage(save_geotiff = FALSE) |> imgRaster()

    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    # Also see if spotPoly is aligned
    v2 <- terra::extract(terra::mean(img), st_centroid(spotPoly(sfe2)))
    expect_true(abs(cor(nCounts, v2$mean)) > 0.4)

    int1 <- st_intersects(rg, spotPoly(sfe))
    int2 <- st_intersects(rowGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)

    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Affine transformation of SFE object with image, after cropping", {
    M <- matrix(c(0.6, -0.2, 0.2, 0.3), nrow = 2)
    v <- c(0, 300)
    sfe <- sfe[,sfe$in_tissue]
    sfe2 <- SpatialFeatureExperiment::affine(sfe, M = M, v = v)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2) |> toSpatRasterImage(save_geotiff = FALSE) |> imgRaster()

    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)

    # Also see if spotPoly is aligned
    v2 <- terra::extract(terra::mean(img), st_centroid(spotPoly(sfe2)))
    expect_true(abs(cor(nCounts, v2$mean)) > 0.4)

    int1 <- st_intersects(rg, spotPoly(sfe))
    int2 <- st_intersects(rowGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)

    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

test_that("Transformation when there's 3D geometry", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))

    sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE,
                     z_option = "3d")
    sfe2 <- SpatialFeatureExperiment::rotate(sfe, degrees = 90)
    ints <- st_intersects(cellSeg(sfe), txSpots(sfe))
    ints2 <- st_intersects(cellSeg(sfe2), txSpots(sfe2))
    expect_equal(ints, ints2)
    rg <- txSpots(sfe2)
    expect_equal(unclass(st_z_range(rg)), c(zmin = 0, zmax = 1),
                 ignore_attr = "crs")
    unlink(dir_use, recursive = TRUE)
})

test_that("Translate SFE object with image", {
    v <- c(100, 200)
    sfe2 <- translate(sfe, v = v)
    bbox_tr <- bbox(sfe2)
    bbox_exp <- bbox(sfe)
    bbox_exp[c("xmin", "xmax")] <- bbox_exp[c("xmin", "xmax")] + v[1]
    bbox_exp[c("ymin", "ymax")] <- bbox_exp[c("ymin", "ymax")] + v[2]
    expect_equal(bbox_tr, bbox_exp)

    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- imgRaster(getImg(sfe2))
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_cg <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(rowGeometry(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, spotPoly(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), spotPoly(sfe2))
    expect_equal(int1, int2)
})

# Affine transform of entire SFE object, for BioFormatsImage=========
library(RBioFormats)
library(EBImage)
fp <- tempdir()
xenium_path <- XeniumOutput(file_path = file.path(fp, "xenium_test"))
try(sfe <- readXenium(xenium_path))
sfe <- readXenium(xenium_path, add_molecules = TRUE)
set.seed(29)
annotGeometry(sfe, "foo") <- ag <-
    st_sf(
        geometry = st_sample(st_bbox(cellSeg(sfe)), 50) |>
            st_buffer(30),
        sample_id = "sample01")

test_that("Transpose SFE object with BioFormatsImage", {
    sfe2 <- SpatialFeatureExperiment::transpose(sfe)
    img <- toExtImage(getImg(sfe2), resolution = 1L)
    # Create image mask
    mask <- img > 500
    spi <- ExtImage(mask, ext = ext(getImg(sfe2))) |> toSpatRasterImage(save_geotiff = FALSE)
    # Due to the way XOA v1 segmentation works, the cell centroid is often
    # outside the nucleus. So I use nuclei centroids
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal(bbox_img_orig[["ymax"]] - bbox_cg_orig[["ymax"]],
                 bbox_cg[["xmin"]] - bbox_img[["xmin"]])

    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with BFI, vertical", {
    sfe2 <- mirror(sfe, direction = "vertical")
    img <- imgRaster(getImg(sfe2), resolution = 1L)
    mask <- img > 500
    spi <- ExtImage(mask, ext = ext(getImg(sfe2))) |> toSpatRasterImage(save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal(bbox_img_orig[["ymax"]] - bbox_cg_orig[["ymax"]],
                 bbox_cg[["ymin"]] - bbox_img[["ymin"]])

    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with BFI, horizontal", {
    sfe2 <- mirror(sfe, direction = "horizontal")
    img <- imgRaster(getImg(sfe2), resolution = 1L)
    mask <- img > 500
    spi <- ExtImage(mask, ext = ext(getImg(sfe2))) |> toSpatRasterImage(save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal(bbox_img_orig[["xmax"]] - bbox_cg_orig[["xmax"]],
                 bbox_cg[["xmin"]] - bbox_img[["xmin"]])

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Rotate SFE object with BFI", {
    sfe2 <- SpatialFeatureExperiment::rotate(sfe, degrees = 45)
    img <- imgRaster(getImg(sfe2), resolution = 1L)
    mask <- img > 500
    spi <- ExtImage(mask, ext = ext(getImg(sfe2))) |> toSpatRasterImage(save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Scale SFE object with BFI", {
    sfe2 <- SpatialFeatureExperiment::scale(sfe, factor = 1.5)
    img <- imgRaster(getImg(sfe2), resolution = 1L)
    mask <- img > 500
    spi <- ExtImage(mask, ext = ext(getImg(sfe2))) |> toSpatRasterImage(save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal((bbox_cg_orig[["xmax"]] - bbox_img_orig[["xmax"]])*1.5,
                 bbox_cg[["xmax"]] - bbox_img[["xmax"]])

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("General affine transformation of SFE object with BFI", {
    M <- matrix(c(0.6, -0.2, 0.2, 0.3), nrow = 2)
    v <- c(0, 300)
    sfe2 <- SpatialFeatureExperiment::affine(sfe, M = M, v = v)
    img <- imgRaster(getImg(sfe2), resolution = 1L)
    mask <- img > 500
    spi <- ExtImage(mask, ext = ext(getImg(sfe2))) |> toSpatRasterImage(save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

# Affine transform of entire SFE object, for ExtImage=========
img <- toExtImage(getImg(sfe), resolution = 1L)
sfe <- addImg(sfe, img, image_id = "exi")

test_that("Transpose SFE object with ExtImage", {
    sfe2 <- SpatialFeatureExperiment::transpose(sfe)
    # Create image mask
    mask <- getImg(sfe2, image_id = "exi") > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    # Due to the way XOA v1 segmentation works, the cell centroid is often
    # outside the nucleus. So I use nuclei centroids
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal(bbox_img_orig[["ymax"]] - bbox_cg_orig[["ymax"]],
                 bbox_cg[["xmin"]] - bbox_img[["xmin"]])

    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with ExtImage, vertical", {
    sfe2 <- mirror(sfe, direction = "vertical")

    mask <- getImg(sfe2, image_id = "exi") > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal(bbox_img_orig[["ymax"]] - bbox_cg_orig[["ymax"]],
                 bbox_cg[["ymin"]] - bbox_img[["ymin"]])

    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Mirror SFE object with ExtImage, horizontal", {
    sfe2 <- mirror(sfe, direction = "horizontal")
    mask <- getImg(sfe2, image_id = "exi") > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal(bbox_img_orig[["xmax"]] - bbox_cg_orig[["xmax"]],
                 bbox_cg[["xmin"]] - bbox_img[["xmin"]])

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Rotate SFE object with ExtImage", {
    sfe2 <- SpatialFeatureExperiment::rotate(sfe, degrees = 45)
    mask <- getImg(sfe2, image_id = "exi") > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("Scale SFE object with ExtImage", {
    sfe2 <- SpatialFeatureExperiment::scale(sfe, factor = 1.5)
    mask <- getImg(sfe2, image_id = "exi") > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_cg_orig <- st_bbox(cellSeg(sfe))
    bbox_img_orig <- ext(getImg(sfe))
    bbox_cg <- st_bbox(cellSeg(sfe2))
    bbox_img <- ext(getImg(sfe2))
    expect_equal((bbox_cg_orig[["xmax"]] - bbox_img_orig[["xmax"]])*1.5,
                 bbox_cg[["xmax"]] - bbox_img[["xmax"]])

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})

test_that("General affine transformation of SFE object with ExtImage", {
    M <- matrix(c(0.6, -0.2, 0.2, 0.3), nrow = 2)
    v <- c(0, 300)
    sfe2 <- SpatialFeatureExperiment::affine(sfe, M = M, v = v)
    mask <- getImg(sfe2, image_id = "exi") > 500
    spi <- toSpatRasterImage(mask, save_geotiff = FALSE)
    v <- terra::extract(spi, st_centroid(nucSeg(sfe2)))
    expect_true(mean(v$lyr.1) > 0.9)

    bbox_rg <- st_bbox(txSpots(sfe2)) |> st_as_sfc()
    bbox_cg <- st_bbox(cellSeg(sfe2)) |> st_as_sfc()
    expect_true(st_covers(bbox_cg, bbox_rg, sparse = FALSE))
    # For ag, can get indices of intersection
    int1 <- st_intersects(ag, cellSeg(sfe))
    int2 <- st_intersects(annotGeometry(sfe2), cellSeg(sfe2))
    expect_equal(int1, int2)
})
# Final cleanup in case failed test messed with cleanup
fp <- tempdir()
unlink(file.path(fp, "cosmx_test"), recursive = TRUE)
unlink(file.path(fp, "vizgen_test"), recursive = TRUE)
unlink(file.path(fp, "xenium_test"), recursive = TRUE)
