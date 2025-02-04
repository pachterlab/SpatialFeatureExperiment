library(S4Vectors)
library(spdep)
set.SubgraphOption(FALSE)
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
    expect_equal(nrow(int_metadata(sfe2)$annotGeometries$baz), 0)
})

test_that("row and col graphs are dropped if drop = TRUE", {
    withr::local_options(SFE_graph_subset = FALSE)
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
    withr::local_options(SFE_graph_subset = FALSE)
    sfe_visium <- sfe_visium[, -1]
    expect_equal(colGraph(sfe_visium, sample_id = "sample01"), g_sub,
        ignore_attr = TRUE
    )
})

test_that("Correctly subset graphs", {
    # Remove one item from sample01
    sfe_visium <- sfe_visium[, -1]
    expect_equal(colGraph(sfe_visium, sample_id = "sample01"), g_sub,
                 ignore_attr = TRUE
    )
})

test_that("Subset the graph when distance-based edge weights are used", {
    attr(g_visium$weights, "mode") <- "distance"
    colGraph(sfe_visium, "foo", "sample01") <- g_visium
    sfe_visium <- sfe_visium[, -1]
    expect_equal(colGraph(sfe_visium, sample_id = "sample01"), g_sub,
                 ignore_attr = TRUE
    )
})

test_that("When only one cell/spot is left in a sample after subsetting", {
    expect_message(sfe1 <- sfe_visium[,1:6], "graphs are meaningless, dropping graphs")
    expect_null(colGraphNames(sfe1, sample_id = "sample02"))
})

test_that("When only one cell is left, option reconstruct", {
    withr::local_options(SFE_graph_subset = FALSE)
    expect_message(sfe1 <- sfe_visium[,1:6], "graphs are meaningless, dropping graphs")
    expect_null(colGraphNames(sfe1, sample_id = "sample02"))
})

test_that("Warning message and dropping graphs when reconstruction info is unavailable", {
    withr::local_options(SFE_graph_subset = FALSE)
    # Remove one item from sample02
    expect_warning(
        sfe_visium <- sfe_visium[, -13],
        "Graph reconstruction info is missing for sample sample02 colGraph bar"
    )
    expect_error(colGraph(sfe_visium, "bar", sample_id = "sample02"))
})

test_that("Warning message and dropping graphs when package required for reconstruction is not installed", {
    withr::local_options(SFE_graph_subset = FALSE)
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

if (!dir.exists("kidney")) dir.create(file.path("kidney", "outs"), recursive = TRUE)
mat_fn <- file.path("kidney", "outs", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_filtered_feature_bc_matrix.h5",
                  destfile = file.path("kidney", "outs", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("kidney", "outs", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_spatial.tar.gz",
                  destfile = file.path("kidney", "outs", "spatial.tar.gz"))
    untar(file.path("kidney", "outs", "spatial.tar.gz"), exdir = file.path("kidney", "outs"))
}

library(sf)
library(SpatialExperiment)
library(terra)
library(SingleCellExperiment)
library(S4Vectors)
sfe1 <- read10xVisiumSFE("ob", sample_id = "ob", unit = "micron")
sfe2 <- read10xVisiumSFE("kidney", sample_id = "kidney", zero.policy = TRUE,
                         unit = "micron")

genes_use <- intersect(rownames(sfe1), rownames(sfe2))
sfe1 <- sfe1[genes_use,]
sfe2 <- sfe2[genes_use,]
rowData(sfe2)$symbol <- rowData(sfe1)$symbol
sfe <- cbind(sfe1, sfe2)

# Cropping: remove the singletons in the kidney dataset
g <- colGraph(sfe2, "visium")
not_singleton <- card(g$neighbour) >= 1
inds <- c(rep(TRUE, ncol(sfe1)), not_singleton)

test_that("Images are cropped after subsetting, multiple samples", {
    sfe3 <- sfe[,inds]
    cg <- st_centroid(st_geometry(spotPoly(sfe3, sample_id = "all")))
    nCounts <- Matrix::colSums(counts(sfe3))
    # For sample 1
    img1 <- getImg(sfe3, sample_id = "ob")
    bbox_geom <- st_bbox(spotPoly(sfe3, "ob")) |> st_as_sfc()
    bbox_img <- ext(img1) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))

    # For sample 2
    img2 <- getImg(sfe3, sample_id = "kidney")
    bbox_geom <- st_bbox(spotPoly(sfe3, "kidney")) |> st_as_sfc()
    bbox_img <- ext(img2) |> st_bbox() |> st_as_sfc()
    expect_true(st_covered_by(bbox_geom, bbox_img, sparse = FALSE))
})

test_that("Images are not cropped when only subsetting rows not cols", {
    sfe3 <- sfe[123,]
    expect_equal(bbox(sfe3, include_image = TRUE), bbox(sfe, include_image = TRUE))
})

sfe1 <- sfe1[rowSums(counts(sfe1)) > 0,]
bbox_use <- bbox(sfe1) |> st_bbox() |> st_as_sfc()
set.seed(29)
ag <- st_sample(bbox_use, 20) |> st_buffer(dist = 100)
ag <- st_sf(geometry = ag, sample_id = "ob")
rg <- st_sample(bbox_use, nrow(sfe1))
rg <- st_sf(geometry = rg)
rownames(rg) <- rownames(sfe1)
annotGeometry(sfe1, "foo") <- ag
rowGeometry(sfe1, "bar") <- rg

test_that("When returning SFE object with 0 rows", {
    # colData shouldn't be affected
    sfe0 <- sfe1[integer(0),]
    expect_equal(dim(sfe0), c(0, ncol(sfe1)))
    expect_equal(names(rowData(sfe0)), names(rowData(sfe1)))
    expect_equal(nrow(rowData(sfe0)), 0)
    expect_equal(rowGeometryNames(sfe0), rowGeometryNames(sfe1))
    expect_equal(colData(sfe0), colData(sfe1))
    expect_equal(annotGeometries(sfe0), annotGeometries(sfe1))
    expect_equal(colGeometries(sfe0), colGeometries(sfe1))
    # Don't want to deal with the SpatRaster pointers
    expect_equal(imgData(sfe0)[,c("sample_id", "image_id", "scaleFactor")],
                 imgData(sfe1)[,c("sample_id", "image_id", "scaleFactor")])
    expect_equal(spatialGraphs(sfe0), spatialGraphs(sfe1))
})

test_that("Returning 0 columns", {
    # Should drop everything with sample_ids
    sfe0 <- sfe1[,integer(0)]
    expect_equal(dim(sfe0), c(nrow(sfe1), 0))
    expect_equal(rowData(sfe0), rowData(sfe1))
    expect_equal(names(colData(sfe0)), names(colData(sfe1)))
    expect_equal(nrow(colData(sfe0)), 0)
    expect_equal(nrow(imgData(sfe0)), 0)
    expect_true(isEmpty(spatialGraphs(sfe0)))
    expect_equal(annotGeometryNames(sfe0), annotGeometryNames(sfe1))
    expect_equal(nrow(annotGeometry(sfe0, "foo")), 0)
})

test_that("Still works when using logical vector of all FALSE", {
    sfe0 <- sfe1[,rep(FALSE, ncol(sfe1))]
    expect_equal(dim(sfe0), c(nrow(sfe1), 0))
    expect_equal(rowData(sfe0), rowData(sfe1))
    expect_equal(names(colData(sfe0)), names(colData(sfe1)))
    expect_equal(nrow(colData(sfe0)), 0)
    expect_equal(nrow(imgData(sfe0)), 0)
    expect_true(isEmpty(spatialGraphs(sfe0)))
    expect_equal(annotGeometryNames(sfe0), annotGeometryNames(sfe1))
    expect_equal(nrow(annotGeometry(sfe0, "foo")), 0)
})
