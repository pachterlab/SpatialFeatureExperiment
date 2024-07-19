library(SFEData)
library(SpatialExperiment)
library(sf)

fp <- tempfile()
fn <- XeniumOutput("v2", file_path = fp)
grid <- st_make_grid(x = st_as_sfc(st_bbox(c(xmin=0, xmax=1027, ymin=4, ymax=1009))),
                     cellsize = 50)
test_that("Directly call aggregateTx to aggregate from file, specify `by`", {
    tx_agged <- aggregateTx(file.path(fn, "transcripts.parquet"), by = grid,
                            spatialCoordsNames = c("x_location", "y_location"),
                            gene_col = "feature_name")
    expect_s4_class(tx_agged, "SpatialFeatureExperiment")
    tx_agged$nCounts <- colSums(counts(tx_agged))
    # empty cells are removed
    expect_true(all(tx_agged$nCounts > 0))
    expect_true(all(st_area(colGeometry(tx_agged)) == 2500))
})

test_that("aggregateTx from file, generate grid", {
    tx_agged <- aggregateTx(file.path(fn, "transcripts.parquet"),
                            spatialCoordsNames = c("x_location", "y_location"),
                            gene_col = "feature_name", cellsize = 50)
    expect_s4_class(tx_agged, "SpatialFeatureExperiment")
    tx_agged$nCounts <- colSums(counts(tx_agged))
    # empty cells are removed
    expect_true(all(tx_agged$nCounts > 0))
    expect_true(all(st_area(colGeometry(tx_agged)) == 2500))
})

test_that("Call aggregateTx for a data frame", {
    df <- read_parquet(file.path(fn, "transcripts.parquet"))
    tx_agged <- aggregateTx(df = df,
                            spatialCoordsNames = c("x_location", "y_location"),
                            gene_col = "feature_name", cellsize = 50)
    expect_s4_class(tx_agged, "SpatialFeatureExperiment")
    tx_agged$nCounts <- colSums(counts(tx_agged))
    # empty cells are removed
    expect_true(all(tx_agged$nCounts > 0))
    expect_true(all(st_area(colGeometry(tx_agged)) == 2500))
})

fn_vizgen <- VizgenOutput("cellpose", file_path = fp)

test_that("aggregateTxTech for Vizgen", {
    sfe <- aggregateTxTech(fn_vizgen, tech = "Vizgen", cellsize = 20)
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), "bins")
    sfe$nCounts <- colSums(counts(sfe))
    # empty cells are removed
    expect_true(all(sfe$nCounts > 0))
    expect_true(all(st_area(colGeometry(sfe)) == 400))
    # Image is aligned
    ids <- imgData(sfe)
    expect_true(nrow(ids) > 1L)
    grid_bbox <- bbox(sfe) |> st_bbox() |> st_as_sfc()
    img_bbox <- bbox(sfe, include_image = TRUE) |> st_bbox() |> st_as_sfc()
    expect_true(st_contains(img_bbox, grid_bbox, sparse = FALSE))
    # I know that this part shouldn't have transcripts
    bbox_no_tx <- c(xmin=6500, xmax = 6539, ymin=-1290, ymax=-1270) |>
        st_bbox() |> st_as_sfc()
    expect_false(any(st_intersects(bbox_no_tx, colGeometry(sfe), sparse = FALSE)))
})

test_that("aggregateTxTech for Xenium", {
    # RBioFormats error
    try(sfe <- aggregateTxTech(fn, tech = "Xenium", cellsize = 50))
    sfe <- aggregateTxTech(fn, tech = "Xenium", cellsize = 50)
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(colGeometryNames(sfe), "bins")
    sfe$nCounts <- colSums(counts(sfe))
    # empty cells are removed
    expect_true(all(sfe$nCounts > 0))
    expect_true(all(st_area(colGeometry(sfe)) == 2500))
    # Image is aligned
    ids <- imgData(sfe)
    expect_true(nrow(ids) > 0)
    grid_bbox <- bbox(sfe) |> st_bbox() |> st_as_sfc()
    img_bbox <- bbox(sfe, include_image = TRUE) |> st_bbox() |> st_as_sfc()
    expect_true(st_contains(img_bbox, grid_bbox, sparse = FALSE))
})

test_that("aggregate for SFE, manually supply `by` argument", {
    sfe <- readXenium(fn)

})

test_that("aggregate for SFE, generate grid from arguments", {

})

test_that("Error message for unacceptable FUN", {

})

test_that("aggregate for SFE, use rowGeometry", {

})

test_that("aggregate for SFE, multiple samples", {

})

unlink(fn, recursive = TRUE)
unlink(fn_vizgen, recursive = TRUE)
