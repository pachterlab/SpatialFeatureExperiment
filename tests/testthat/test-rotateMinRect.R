library(sf)

g <- structure(list(structure(c(4399.55012508914, 3803.07118317876
), class = c("XY", "POINT", "sfg")), structure(c(2431.08190332779, 
2292.70498132735), class = c("XY", "POINT", "sfg")), structure(c(3347.88654099145, 
3266.9065113657), class = c("XY", "POINT", "sfg")), structure(c(3707.26120687946, 
3663.17437017167), class = c("XY", "POINT", "sfg")), structure(c(3750.39171528186, 
4017.8412437251), class = c("XY", "POINT", "sfg")), structure(c(776.023716388489, 
888.486086093483), class = c("XY", "POINT", "sfg")), structure(c(3589.0257674615, 
3471.39714713708), class = c("XY", "POINT", "sfg")), structure(c(1125.6967871746, 
941.467827825188), class = c("XY", "POINT", "sfg")), structure(c(4195.76846125686, 
3801.62926699026), class = c("XY", "POINT", "sfg")), structure(c(3240.8615299643, 
2742.67244338588), class = c("XY", "POINT", "sfg"))), class = c("sfc_POINT", 
"sfc"), precision = 0, bbox = structure(c(xmin = 776.023716388489, 
ymin = 888.486086093483, xmax = 4399.55012508914, ymax = 4017.8412437251
), class = "bbox"), crs = structure(list(input = NA_character_, 
    wkt = NA_character_), class = "crs"), n_empty = 0L)

test_that("rotateMinRect for matrix", {
    mat <- st_coordinates(g)[,1:2]
    r <- rotateMinRect(mat, orientation = "horizontal")
    expect_equal(class(r), c("matrix", "array"))
    r <- st_multipoint(r)
    rc <- st_minimum_rotated_rectangle(st_union(r)) |> st_coordinates()
    bb <- st_bbox(r)
    expect_true(bb["xmax"]-bb["xmin"] >  bb["ymax"]-bb["ymin"])
    # That it's not at away from 0 or 90 degrees
    expect_true(anyDuplicated(rc[,1])>0)
})

test_that("rotateMinRect for sfc, horizontal", {
    r <- rotateMinRect(g, orientation = "horizontal")
    expect_s3_class(r, "sfc")
    rc <- st_minimum_rotated_rectangle(st_union(r)) |> st_coordinates()
    bb <- st_bbox(r)
    expect_true(bb["xmax"]-bb["xmin"] >  bb["ymax"]-bb["ymin"])
    # That it's not at away from 0 or 90 degrees
    expect_true(anyDuplicated(rc[,1])>0)
    # Centroid didn't shift a lot
    c1 <- st_bbox(g) |> st_as_sfc()
    c2 <- bb |> st_as_sfc()
    expect_true(st_intersects(c1, c2, sparse = FALSE))
})

test_that("rotateMinRect for sfc, vertical", {
    r <- rotateMinRect(g, orientation = "vertical")
    rc <- st_minimum_rotated_rectangle(st_union(r)) |> st_coordinates()
    bb <- st_bbox(r)
    expect_true(bb["xmax"]-bb["xmin"] <  bb["ymax"]-bb["ymin"])
    # That it's not at away from 0 or 90 degrees
    expect_true(anyDuplicated(rc[,1])>0)
    # Centroid didn't shift a lot
    c1 <- st_bbox(g) |> st_as_sfc()
    c2 <- bb |> st_as_sfc()
    expect_true(st_intersects(c1, c2, sparse = FALSE))
})

test_that("rotateMinRect for sf", {
    df <- st_sf(geometry = g, crs = NA)
    r <- rotateMinRect(df, orientation = "horizontal")
    expect_s3_class(r, "sf")
    rc <- st_minimum_rotated_rectangle(st_union(r)) |> st_coordinates()
    bb <- st_bbox(r)
    expect_true(bb["xmax"]-bb["xmin"] >  bb["ymax"]-bb["ymin"])
    # That it's not at away from 0 or 90 degrees
    expect_true(anyDuplicated(rc[,1])>0)
})

test_that("rotateMinRect for SFE", {
    sfe <- SpatialFeatureExperiment(list(foo = matrix(rnorm(20), nrow = 2, ncol = 10)),
                                    colGeometries = list(centroids = st_sf(geometry = g)))
    colnames(sfe) <- 1:10
    sfe2 <- rotateMinRect(sfe, "horizontal")
    rc <- st_minimum_rotated_rectangle(st_union(centroids(sfe2))) |> st_coordinates()
    bb <- st_bbox(centroids(sfe2))
    expect_true(bb["xmax"]-bb["xmin"] >  bb["ymax"]-bb["ymin"])
    # That it's not at away from 0 or 90 degrees
    expect_true(anyDuplicated(rc[,1])>0)
})

