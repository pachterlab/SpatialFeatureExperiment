library(sf)

sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
    package = "SpatialFeatureExperiment"
))
sfe_visium <- addVisiumSpotPoly(sfe_visium, 0.5)
bbox_use <- st_as_sfc(st_bbox(c(xmin = 1, xmax = 3, ymin = 1, ymax = 6),
    crs = NA
))
bbox_use2 <- st_sf(geometry = bbox_use, sample_id = "sample01", crs = NA)

test_that("All spots in the cropped SFE objects indeed are covered by the bbox", {
    sfe_cropped <- crop(sfe_visium, bbox_use, sample_id = "all")
    cg <- spotPoly(sfe_cropped, "all")
    expect_true(all(st_any_pred(cg, bbox_use, pred = st_covered_by)))
    expect_true(st_geometry_type(cg, by_geometry = FALSE) == "POLYGON")
    sfe_cropped2 <- crop(sfe_visium, y = bbox_use2, sample_id = "sample01")
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

test_that("When a geometry is broken into multiple pieces", {
    notch <- st_as_sfc(st_bbox(c(xmin = 0, xmax = 1.5, ymin = 1.7, ymax = 1.9)))
    bbox_use3 <- st_difference(bbox_use, notch)
    sfe_cropped3 <- crop(sfe_visium, bbox_use3, sample_id = "sample01")
    cg <- spotPoly(sfe_cropped3, "all")
    expect_true(st_geometry_type(cg, by_geometry = FALSE) == "MULTIPOLYGON")
})

annotGeometry(sfe_visium, "bbox", sample_id = "sample01") <- bbox_use2
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

test_that("Find bbox of samples", {
    cg1 <- spotPoly(sfe_visium, "sample01")
    bbox1 <- bbox(sfe_visium, "sample01")
    expect_equal(st_bbox(st_union(cg1, bbox_use)), bbox1,
        ignore_attr = c("class", "crs")
    )
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

ag1 <- st_sfc(st_linestring(matrix(c(1, 5, 2, 2), ncol = 2, byrow = TRUE)))
ag2 <- st_sfc(st_linestring(matrix(c(2,2,1,5), ncol = 2, byrow = TRUE)))
ag <- c(ag1, ag2)
ag <- st_sf(sample_id = c("sample01", "sample02"),
            geometry = ag)
annotGeometry(sfe_visium, "foo", sample_id = "all") <- ag
test_that("Transpose SFE, no images", {
    sfe_visium2 <- transpose(sfe_visium)
    coords1 <- spatialCoords(sfe_visium)
    coords2 <- spatialCoords(sfe_visium2)
    expect_equal(coords1[,1], coords2[,2])
    coords_cg1 <- st_coordinates(spotPoly(sfe_visium, sample_id = "all"))
    coords_cg2 <- st_coordinates(spotPoly(sfe_visium2, sample_id = "all"))
    expect_equal(coords_cg1[,1], coords_cg2[,2])
    coords_ag1 <- st_coordinates(annotGeometry(sfe_visium, "foo", "all"))
    coords_ag2 <- st_coordinates(annotGeometry(sfe_visium2, "foo", "all"))
    expect_equal(coords_ag1[,1], coords_ag2[,2])
})

test_that("Mirror SFE, no images, vertical", {
    sfe_visium2 <- mirror(sfe_visium)
    coords1 <- spatialCoords(sfe_visium)
    coords2 <- spatialCoords(sfe_visium2)
    expect_equal(coords1[,2], -coords2[,2])
    expect_equal(coords1[,1], coords2[,1])
    coords_cg1 <- st_coordinates(spotPoly(sfe_visium, sample_id = "all"))
    coords_cg2 <- st_coordinates(spotPoly(sfe_visium2, sample_id = "all"))
    expect_equal(coords_cg1[,2], -coords_cg2[,2])
    expect_equal(coords_cg1[,1], coords_cg2[,1])
    coords_ag1 <- st_coordinates(annotGeometry(sfe_visium, "foo", "all"))
    coords_ag2 <- st_coordinates(annotGeometry(sfe_visium2, "foo", "all"))
    expect_equal(coords_ag1[,1], coords_ag2[,1])
    expect_equal(coords_ag1[,2], -coords_ag2[,2])
})

test_that("Mirror SFE, no images, horizontal", {
    sfe_visium2 <- mirror(sfe_visium, direction = "horizontal")
    coords1 <- spatialCoords(sfe_visium)
    coords2 <- spatialCoords(sfe_visium2)
    expect_equal(coords1[,1], -coords2[,1])
    expect_equal(coords1[,2], coords2[,2])
    coords_cg1 <- st_coordinates(spotPoly(sfe_visium, sample_id = "all"))
    coords_cg2 <- st_coordinates(spotPoly(sfe_visium2, sample_id = "all"))
    expect_equal(coords_cg1[,2], coords_cg2[,2])
    expect_equal(coords_cg1[,1], -coords_cg2[,1])
    coords_ag1 <- st_coordinates(annotGeometry(sfe_visium, "foo", "all"))
    coords_ag2 <- st_coordinates(annotGeometry(sfe_visium2, "foo", "all"))
    expect_equal(coords_ag1[,1], -coords_ag2[,1])
    expect_equal(coords_ag1[,2], coords_ag2[,2])
})

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

sfe_shifted <- removeEmptySpace(sfe)
bbox_new <- st_as_sfc(st_bbox(bbox(sfe_shifted, sample_id = "Vis5A")))
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

test_that("annotSummary", {
    out <- annotSummary(sfe1, "spotPoly", "myofiber_simplified", "area")
    expect_s3_class(out, "data.frame")
    expect_equal(rownames(out), colnames(sfe1))
    expect_equal(names(out), "area")
    expect_true(is.numeric(out$area))
})

# Need uncropped image
if (!dir.exists("outs")) dir.create("outs")
mat_fn <- file.path("outs", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_filtered_feature_bc_matrix.h5",
                  destfile = file.path("outs", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("outs", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_Mouse_Olfactory_Bulb/Visium_Mouse_Olfactory_Bulb_spatial.tar.gz",
                  destfile = file.path("outs", "spatial.tar.gz"))
    untar(file.path("outs", "spatial.tar.gz"), exdir = "outs")
}
library(sf)
library(SpatialExperiment)
library(terra)
library(SingleCellExperiment)
sfe <- read10xVisiumSFE(".")
test_that("Image is shifted after removing empty space", {
    sfe2 <- removeEmptySpace(sfe)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
    expect_true(st_area(bbox_geom) / st_area(bbox_img) > 0.98)
})

test_that("Image is cropped after cropping SFE object", {
    bc <- bbox(sfe)
    bbox_use <- c(xmin = bc["xmin"], xmax = bc["xmin"] + 2000,
                  ymin = bc["ymin"], ymax = bc["ymin"] + 2000) |>
        setNames(c("xmin", "xmax", "ymin", "ymax")) |>
        st_bbox() |> st_as_sfc()
    sfe2 <- SpatialFeatureExperiment::crop(sfe, bbox_use)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean, use = "complete.obs")) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
    expect_true(st_area(bbox_geom) / st_area(bbox_img) > 0.98)
})

test_that("Transpose SFE object with image", {
    sfe2 <- transpose(sfe)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
})

test_that("Transpose SFE object with image, after cropping image", {
    sfe <- sfe[,sfe$in_tissue]
    sfe2 <- transpose(sfe)
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
})

test_that("Mirror SFE object with image, vertical", {
    sfe2 <- mirror(sfe, direction = "vertical")
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
})

test_that("Mirror SFE object with image, horizontal", {
    sfe2 <- mirror(sfe, direction = "horizontal")
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
})
