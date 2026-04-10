library(SFEData)
library(sf)

# Tissue boundary from brightfield image=======================
if (!dir.exists("kidney")) dir.create("kidney", recursive = TRUE)
mat_fn <- file.path("kidney", "filtered_feature_bc_matrix.h5")
if (!file.exists(mat_fn))
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_filtered_feature_bc_matrix.h5",
                  destfile = file.path("kidney", "filtered_feature_bc_matrix.h5"),
                  mode = "wb")
if (!dir.exists(file.path("kidney", "spatial"))) {
    download.file("https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_spatial.tar.gz",
                  destfile = file.path("kidney", "spatial.tar.gz"))
    untar(file.path("kidney", "spatial.tar.gz"), exdir = "kidney")
}
sfe <- read10xVisiumSFE(dirs = "kidney", zero.policy = TRUE)

test_that("Find tissue boundary from H&E", {
    tb <- getTissueBoundaryImg(sfe, image_id = "hires", fill_holes = TRUE,
                               image_type = "brightfield")
    expect_s3_class(tb, "sf")
    expect_equal(names(tb), c("sample_id", "geometry"))
    expect_equal(tb$sample_id, "sample01")
    expect_equal(nrow(tb), 1L)
    prop_in_tissue <- mean(st_intersects(tb, spotPoly(sfe), sparse = FALSE))
    expect_true(prop_in_tissue > 0.95) # there's debris so not 100%
})

# Tissue boundary from fluorescent image=======================
# Note that it doesn't work well when the image looks sparse
fp <- tempfile()
fn <- VizgenOutput("cellpose", file_path = fp)
suppressWarnings(sfe <- readVizgen(fn))

test_that("Find tissue boundary from fluorescent image", {
    tb <- getTissueBoundaryImg(sfe, image_id = "PolyT_z2",
                               image_type = "fluorescent", n_pieces = 1)
    expect_equal(names(tb), c("sample_id", "geometry"))
    expect_equal(tb$sample_id, "sample01")
    expect_equal(nrow(tb), 1L)
    prop_in_tissue <- mean(st_intersects(tb, cellSeg(sfe), sparse = FALSE))
    expect_true(prop_in_tissue > 0.95)
})

# Tissue boundary from concave hull================
test_that("Concave hull tissue boundary, 1 piece", {
    tb <- getTissueBoundaryConcave(sfe, ratio = 0.12)
    expect_equal(names(tb), c("sample_id", "geometry"))
    expect_equal(tb$sample_id, "sample01")
    expect_equal(nrow(tb), 1L)
    prop_in_tissue <- mean(st_intersects(tb, cellSeg(sfe), sparse = FALSE))
    expect_true(prop_in_tissue > 0.95)
})

fn <- XeniumOutput("v2", file_path = fp)
sfe <- readXenium(fn)
test_that("Concave hull tissue boundary, multiple pieces", {
    tb <- getTissueBoundaryConcave(sfe, multiple_pieces = TRUE,
                                   distance_cutoff = 25, ratio = 0.1)
    expect_equal(names(tb), c("sample_id", "geometry"))
    expect_equal(tb$sample_id, "sample01")
    expect_equal(nrow(tb), 1L)
    expect_equal(st_geometry_type(tb$geometry) |> as.character(), "MULTIPOLYGON")
    g <- st_cast(st_geometry(tb$geometry), "POLYGON")
    expect_equal(length(g), 2L)
    prop_in_tissue <- mean(st_intersects(tb, cellSeg(sfe), sparse = FALSE))
    expect_true(prop_in_tissue > 0.95)
})
