# Unit test the subsetting method
sfe2 <- readRDS(system.file("extdata/sfe_multi_sample.rds",
    package = "SpatialFeatureExperiment"
))
agr1 <- readRDS(system.file("extdata/annotgraph1.rds",
    package = "SpatialFeatureExperiment"
))
ag <- readRDS(system.file("extdata/ag.rds",
    package = "SpatialFeatureExperiment"
))
# Should have passed unit test for graph_wrapper
spatialGraph(sfe2, type = "foo", MARGIN = 2, sample_id = "sample01") <-
    findSpatialNeighbors(sfe2, "sample01",
        type = "spatialCoords", MARGIN = 2,
        method = "tri2nb"
    )
annotGraph(sfe2, "bar", "sample01") <- agr1
annotGeometry(sfe2, "baz", "sample01") <- ag

test_that("After removing one sample_id, it's also removed in annotGeometries", {
    sfe2 <- sfe2[, 4:5]
    expect_true(is.null(int_metadata(sfe2)$annotGeometries$baz))
})

test_that("row and col graphs are dropped if drop = TRUE", {
    expect_message(sfe2 <- sfe2[, 2:5, drop = TRUE], "Dropping all")
    # Don't have rowGraphs to begin with
    expect_true(is.null(unlist(as.list(colGraphs(sfe2)))))
})

sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
    package = "SpatialFeatureExperiment"
))
g_visium <- readRDS(system.file("extdata/colgraph_visium.rds",
    package = "SpatialFeatureExperiment"
))
g_visium2 <- readRDS(system.file("extdata/colgraph_visium2.rds",
    package = "SpatialFeatureExperiment"
))
g_sub <- readRDS(system.file("extdata/colgraph_visium_sub.rds",
    package = "SpatialFeatureExperiment"
))
colGraph(sfe_visium, "foo", "sample01") <- g_visium
colGraph(sfe_visium, "bar", "sample02") <- g_visium2

test_that("Retain correct spatialGraphs structure when one entire sample is left", {
    sfe_visium1 <- sfe_visium[, colData(sfe_visium)$sample_id == "sample01"]
    sgr <- int_metadata(sfe_visium1)[["spatialGraphs"]]
    expect_true(is(sgr, "DataFrame"))
    expect_equal(names(sgr), "sample01")
    expect_true(is(spatialGraph(sfe_visium1, "foo", 2, "sample01"), "listw"))
})

test_that("Correctly reconstruct the graphs when they need to be reconstructed", {
    # Remove one item from sample01
    sfe_visium <- sfe_visium[, -1]
    expect_equal(colGraph(sfe_visium, sample_id = "sample01"), g_sub,
        ignore_attr = c("call", "method")
    )
})

test_that("Warning message and dropping graphs when reconstruction info is unavailable", {
    # Remove one item from sample02
    expect_warning(
        sfe_visium <- sfe_visium[, -13],
        "Graph reconstruction info is missing for sample sample02 colGraph bar"
    )
    expect_error(colGraph(sfe_visium, "bar", sample_id = "sample02"))
})

test_that("Warning message and dropping graphs when package required for reconstruction is not installed", {
    attr(g_visium2, "method") <- list(
        FUN = "findVisiumGraph",
        package = "foobar",
        args = list(
            style = "W",
            zero.policy = NULL,
            sample_id = "sample01"
        )
    )
    colGraph(sfe_visium, "bar", "sample02") <- g_visium2
    expect_warning(
        sfe_visium <- sfe_visium[, -13],
        "Package foobar used to construct graph for sample sample02 colGraph bar is not installed"
    )
    expect_error(colGraph(sfe_visium, "bar", sample_id = "sample02"))
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
test_that("Images are cropped after subsetting", {
    bc <- bbox(sfe)
    bbox_use <- c(xmin = bc["xmin"], xmax = bc["xmin"] + 2000,
                  ymin = bc["ymin"], ymax = bc["ymin"] + 2000) |>
        setNames(c("xmin", "xmax", "ymin", "ymax")) |>
        st_bbox() |> st_as_sfc()
    sfe2 <- sfe[,st_intersects(spotPoly(sfe), bbox_use, sparse = FALSE)[,1]]
    cg <- df2sf(spatialCoords(sfe2), spatialCoordsNames(sfe2))
    img <- getImg(sfe2)@image
    v <- terra::extract(terra::mean(img), cg)
    nCounts <- Matrix::colSums(counts(sfe2))
    expect_true(abs(cor(nCounts, v$mean)) > 0.4)
    bbox_geom <- st_bbox(spotPoly(sfe2)) |> st_as_sfc()
    bbox_img <- as.vector(ext(img)) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
})