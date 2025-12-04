library(SFEData)
library(SpatialExperiment)
library(sparseMatrixStats)
library(sf)
library(arrow)

fp <- tempfile()
fn <- XeniumOutput("v2", file_path = fp)
grid <- st_make_grid(x = st_as_sfc(st_bbox(c(xmin=0, xmax=1027, ymin=4, ymax=1009))),
                     cellsize = 50)
gene_names <- read_parquet(file.path(fn, "transcripts.parquet"), col_select = "feature_name")[[1]] |> unique()
test_that("Directly call aggregateTx to aggregate from file, specify `by`", {
    tx_agged <- aggregateTx(file.path(fn, "transcripts.parquet"), by = grid,
                            spatialCoordsNames = c("x_location", "y_location"),
                            gene_col = "feature_name")
    expect_s4_class(tx_agged, "SpatialFeatureExperiment")
    expect_true(inherits(counts(tx_agged), "matrix"))
    tx_agged$nCounts <- colSums(counts(tx_agged))
    # empty cells are removed
    expect_true(all(tx_agged$nCounts > 0))
    expect_true(all(st_area(colGeometry(tx_agged)) == 2500))
    expect_true(all(rownames(tx_agged) %in% gene_names))
    expect_equal(colnames(tx_agged), as.character(seq_along(grid)))
})

test_that("Directly call aggregateTx, sparse = TRUE", {
    tx_agged <- aggregateTx(file.path(fn, "transcripts.parquet"), by = grid,
                            spatialCoordsNames = c("x_location", "y_location"),
                            gene_col = "feature_name", sparse = TRUE)
    expect_s4_class(tx_agged, "SpatialFeatureExperiment")
    expect_s4_class(counts(tx_agged), "dgCMatrix")
    tx_agged$nCounts <- colSums(counts(tx_agged))
    # empty cells are removed
    expect_true(all(tx_agged$nCounts > 0))
    expect_true(all(rownames(tx_agged) %in% gene_names))
    expect_equal(colnames(tx_agged), as.character(seq_along(grid)))
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
    expect_true(all(rownames(tx_agged) %in% gene_names))
    expect_equal(colnames(tx_agged), as.character(seq_along(grid)))
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
    expect_true(all(rownames(tx_agged) %in% gene_names))
    expect_equal(colnames(tx_agged), as.character(seq_along(grid)))
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

sfe <- readXenium(fn, add_molecules = TRUE)
# Deal with logical and categorical variables in colData
set.seed(29)
sfe$logical <- sample(c(TRUE, FALSE), ncol(sfe), replace = TRUE)
sfe$categorical <- sample(LETTERS[1:3], ncol(sfe), replace = TRUE)

pieces <- readRDS(system.file("extdata/pieces.rds", package = "SpatialFeatureExperiment"))
sfes <- splitByCol(sfe, pieces)
sfes[[2]] <- changeSampleIDs(sfes[[2]], c(sample01 = "sample02"))
sfe2 <- cbind(sfes[[1]], sfes[[2]])

grid2 <- st_make_grid(cellSeg(sfe), cellsize = 50)

test_that("Error messages when checking input to aggregate.SFE", {
    expect_error(aggregate(sfe), "Either `by` or `cellsize` must be specified.")
    expect_error(aggregate(sfe, by = list(foo = grid)), "None of the geometries in `by` correspond to sample_id")
    expect_error(aggregate(sfe, by = iris$Sepal.Length), "`by` must be either sf or sfc")
    expect_error(aggregate(sfe2, by = grid), "`by` must be an sf data frame with a column `sample_id`")
    expect_error(aggregate(sfe, by = grid), "`by` does not overlap with this sample")
})

test_that("aggregate for SFE by cells, manually supply `by` argument", {
    agg <- aggregate(sfe, by = grid2)
    expect_s4_class(agg, "SpatialFeatureExperiment")
    # empty grid cells were removed
    expect_true(ncol(agg) <= length(grid2) & ncol(agg) > 0)
    expect_true(all(colSums(counts(agg)) > 0))
    expect_equal(rowGeometry(sfe), rowGeometry(agg))
    expect_equal(imgData(sfe), imgData(agg))
    expect_type(agg$categorical, "list")
    expect_true(is.numeric(agg$logical))
    expect_true(all(agg$logical >= 0L))
})

test_that("aggregate for SFE by cells, sf `by` argument, single sample", {
    grid_sf <- st_sf(geometry = grid2)
    agg <- aggregate(sfe, by = grid_sf)
    expect_true(ncol(agg) <= length(grid2) & ncol(agg) > 0)
    expect_true(all(colSums(counts(agg)) > 0))
})

test_that("aggregate with sf `by` argument, multiple samples", {
    grid3 <- rbind(st_sf(geometry = grid2, sample_id = "sample01"),
                   st_sf(geometry = grid2, sample_id = "sample02"))
    agg <- aggregate(sfe2, by = grid3)
    expect_true(ncol(agg) <= nrow(grid3) & ncol(agg) > 0)
    expect_true(all(colSums(counts(agg)) > 0))
})

test_that("aggregate when some cells are outside bins", {
    grid4 <- grid2[seq_len(floor(length(grid2)/2))]
    agg <- aggregate(sfe, by = grid4)
    expect_true(ncol(agg) <= length(grid4) & ncol(agg) > 0)
    expect_true(all(colSums(counts(agg)) > 0))
})

test_that("aggregate.SFE use a row* function", {
    agg2 <- aggregate(sfe, by = grid2, FUN = rowMedians)
    expect_s4_class(agg2, "SpatialFeatureExperiment")
    # empty grid cells were removed
    expect_true(ncol(agg2) <= length(grid2) & ncol(agg2) > 0)
    expect_true(all(colSums(counts(agg2)) >= 0))
    expect_equal(rowGeometry(sfe), rowGeometry(agg2))
    expect_equal(imgData(sfe), imgData(agg2))
    expect_type(agg2$categorical, "list")
    expect_true(is.numeric(agg2$logical))
    expect_true(all(agg2$logical >= 0L))
})

test_that("aggregate for SFE, generate grid from arguments", {
    agg <- aggregate(sfe, cellsize = 50)
    expect_s4_class(agg, "SpatialFeatureExperiment")
    # empty grid cells were removed
    expect_true(ncol(agg) <= length(grid2) & ncol(agg) > 0)
    expect_true(all(colSums(counts(agg)) > 0))
    expect_equal(rowGeometry(sfe), rowGeometry(agg))
    expect_equal(imgData(sfe), imgData(agg))
    expect_type(agg$categorical, "list")
    expect_true(is.numeric(agg$logical))
    expect_true(all(agg$logical >= 0L))
})

test_that("aggregate for SFE, use rowGeometry", {
    agg <- aggregate(sfe, by = grid2, rowGeometryName = "txSpots")
    expect_s4_class(agg, "SpatialFeatureExperiment")
    # empty grid cells were removed
    expect_true(ncol(agg) <= length(grid2) & ncol(agg) > 0)
    expect_true(all(colSums(counts(agg)) > 0))
    expect_equal(rowGeometry(sfe), rowGeometry(agg))
    expect_equal(imgData(sfe), imgData(agg))
})

test_that("aggregate for SFE, multiple samples", {
    agg2 <- aggregate(sfe2, cellsize = 50)
    expect_s4_class(agg2, "SpatialFeatureExperiment")
    # empty grid cells were removed
    expect_true(ncol(agg2) <= length(grid2) & ncol(agg2) > 0)
    expect_true(all(colSums(counts(agg2)) >= 0))
    expect_equal(rowGeometries(sfe2), rowGeometries(agg2))
    expect_equal(imgData(sfe2), imgData(agg2))
    expect_type(agg2$categorical, "list")
    expect_true(is.numeric(agg2$logical))
    expect_true(all(agg2$logical >= 0L))
    expect_equal(sampleIDs(sfe2), sampleIDs(agg2))
})

test_that("aggregate for SFE using txSpots, multiple samples", {
    agg2 <- aggregate(sfe2, cellsize = 50, rowGeometryName = "txSpots")
    expect_equal(sampleIDs(agg2), sampleIDs(sfe2))
    expect_equal(nrow(agg2), nrow(sfe2))
    expect_equal(rowGeometries(agg2), rowGeometries(sfe2))
})

test_that("Empty geometries in rowGeometry", {
    sfe <- readXenium(fn, add_molecules = TRUE)
    bbox2 <- c(xmin=700, xmax=850, ymin=-800, ymax=-650)
    sfe_sub <- crop(sfe, bbox2)
    
    txSpots(sfe_sub) <- st_cast(txSpots(sfe_sub), "MULTIPOINT") |> st_zm()
    sfe_agg_cell <- aggregate(sfe_sub, cellSeg(sfe_sub), rowGeometryName = "txSpots")
    expect_s4_class(counts(sfe_agg_cell), "dgCMatrix")
    expect_true(ncol(sfe_agg_cell) <= ncol(sfe_sub)) # some cells don't have transcripts here
    expect_equal(nrow(sfe_agg_cell), nrow(sfe_sub))
    expect_equal(rowGeometries(sfe_agg_cell), rowGeometries(sfe_sub))
})

unlink(fn, recursive = TRUE)
unlink(fn_vizgen, recursive = TRUE)
