library(SFEData)
library(sf)
test_that("Read MERFISH transcript spots into rowGeometries", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    sfe <- readVizgen(dir_use, z = 3L, image = "PolyT", add_molecules = TRUE)
    expect_equal(rowGeometryNames(sfe), "txSpots")
    rg <- txSpots(sfe)
    expect_equal(as.character(st_geometry_type(rg, FALSE)), "MULTIPOINT")
    # Check that the spots are flipped and aligned with the image
    img <- getImg(sfe)
    v <- terra::extract(img, rg)
    expect_true(sum(v$mosaic_PolyT_z3 < 30, na.rm = TRUE) < 10)
    unlink(dir_use, recursive = TRUE)
})

test_that("Format MERFISH transcript spots for colGeometries", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    expect_error(formatTxSpots(file.path(dir_use, "detected_transcripts.csv"),
                               dest = "colGeometry"),
                 "file_out must be specified")
    expect_error(formatTxSpots(file.path(dir_use, "detected_transcripts.csv"),
                               dest = "colGeometry", file_out = "vizgen_cellbound/tx_in_cells"),
                 "Column indicating cell ID not found.")

    df <- read.csv(file.path(dir_use, "detected_transcripts.csv"),
                   header = TRUE)
    sfe <- readVizgen(dir_use, z = 3L, image = "PolyT")
    df$cell_id <- sample(colnames(sfe), nrow(df), replace = TRUE)
    rownames(df) <- df$X
    df$X <- NULL
    write.csv(df, file = file.path(dir_use, "detected_transcripts.csv"),
              row.names = TRUE, quote = FALSE)

    # First run
    cg <- formatTxSpots(file.path(dir_use, "detected_transcripts.csv"),
                        dest = "colGeometry",
                        file_out = file.path(dir_use, "tx_in_cells"),
                        z = 3L)
    dir_check <- file.path(dir_use, "tx_in_cells")
    expect_equal(cg, dir_check)
    expect_true(dir.exists(dir_check))
    fns_expect <- paste0(unique(df$gene[df$global_z == 3L]), "_spots.parquet")
    fns <- list.files(dir_check)
    expect_setequal(fns, fns_expect)

    # Check contents
    fn <- file.path(dir_check, fns_expect[1])
    g <- sfarrow::st_read_parquet(fn)
    expect_equal(as.character(st_geometry_type(g, FALSE)), "MULTIPOINT")

    # Second run
    time_note <- Sys.time() # Check the files weren't written again
    cg <- formatTxSpots(file.path(dir_use, "detected_transcripts.csv"),
                        dest = "colGeometry",
                        file_out = file.path(dir_use, "tx_in_cells"))
    expect_equal(cg, normalizePath(dir_check))
    time_check <- file.info(fn)$mtime
    expect_true(time_note > time_check)

    unlink(dir_use, recursive = TRUE)
})

test_that("Error messages in formatTxSpots", {
    fp <- tempdir()
    dir_use <- VizgenOutput("hdf5", file_path = file.path(fp, "vizgen_test"))
    file <- file.path(dir_use, "detected_transcripts.csv")
    expect_error(formatTxSpots(file, dest = "colGeometry", file_out = NULL),
                 "file_out must be specified")
    file2 <- system.file("extdata/sfe_visium.rds", package = "SpatialFeatureExperiment")
    expect_error(formatTxSpots(file2),
                 "The file must be one of csv")
    expect_error(formatTxSpots(file, z = "foo"),
                 "z must either be numeric")
    expect_error(formatTxSpots(file, spatialCoordsNames = c("foo", "bar")),
                 "foo, bar not found")
    expect_error(formatTxSpots(file, z = 8L),
                 "z plane\\(s\\) specified not found")
    expect_error(formatTxSpots(file, dest = "colGeometry", cell_col = "foo",
                               file_out = "bar"),
                 "Column indicating cell ID not found")
    unlink(dir_use, recursive = TRUE)
})

test_that("Format CosMX spots for colGeometry, multiple z-planes", {
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))

    cg <- formatTxSpots(file.path(dir_use, "Run5642_S3_Quarter_tx_file.csv"),
                        dest = "colGeometry", z = "all", z_option = "split",
                        cell_col = c("cell_ID", "fov"),
                        gene_col = "target", not_in_cell_id = "0",
                        spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
                        file_out = file.path(dir_use, "tx_spots"))
    expect_equal(cg, file.path(dir_use, "tx_spots"))
    # Oh, great, there's Bex1/2, illegal file name. No wonder people don't like CosMX
    df <- data.table::fread(file.path(dir_use, "Run5642_S3_Quarter_exprMat_file.csv"))
    genes <- names(df)[-c(1:2)]
    combs <- expand.grid(gene = genes, z = 0:1, stringsAsFactors = FALSE)
    fns_expect <- paste0(gsub("/", ".", combs$gene), "_spots_z", combs$z, ".parquet")
    fns <- list.files(file.path(dir_use, "tx_spots"))
    # Not all genes have spots in this downsampled toy dataset
    expect_true(all(fns %in% fns_expect))

    fn <- file.path(dir_use, "tx_spots", fns[1])
    g <- sfarrow::st_read_parquet(fn)
    expect_equal(st_geometry_type(g, FALSE) |> as.character(), "MULTIPOINT")
    unlink(dir_use, recursive = TRUE)
})

test_that("Read transcript spots from a subset of genes", {
    skip_if_not(gdalParquetAvailable())
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    fn_tx <- formatTxTech(fn, tech = "Xenium", flip = TRUE, return = FALSE,
                          file_out = file.path(fn, "tx_spots.parquet"))
    gene_select <- c("ACE2", "BMX")
    df <- readSelectTx(fn_tx, gene_select)
    expect_equal(as.character(st_geometry_type(df, FALSE)), "MULTIPOINT")
    expect_equal(nrow(df), 2)
    expect_equal(rownames(df), gene_select)
    expect_equal(names(df), c("gene", "codeword_index", "geometry"))
    unlink(fn, recursive = TRUE)
})

test_that("Error message when Parquet driver is unavailable", {
    skip_if(gdalParquetAvailable())
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    fn_tx <- formatTxTech(fn, tech = "Xenium", flip = TRUE, return = FALSE,
                          file_out = file.path(fn, "tx_spots.parquet"))
    gene_select <- c("ACE2", "BMX")
    expect_error(df <- readSelectTx(fn_tx, gene_select),
                 "GDAL Parquet driver is required")
    unlink(fn, recursive = TRUE)
})

test_that("Add a subset of spots", {
    skip_if_not(gdalParquetAvailable())
    library(RBioFormats)
    fp <- tempdir()
    fn <- XeniumOutput("v2", file_path = file.path(fp, "xenium_test"))
    try(sfe <- readXenium(fn))
    sfe <- readXenium(fn)
    fn_tx <- formatTxTech(fn, tech = "Xenium", flip = TRUE, return = FALSE,
                          file_out = file.path(fn, "tx_spots.parquet"))
    sfe <- addSelectTx(sfe, fn_tx, head(rownames(sfe), 5), swap_rownames = "Symbol")
    expect_equal(as.character(st_geometry_type(txSpots(sfe), FALSE)), "MULTIPOINT")
    is_empty <- st_is_empty(txSpots(sfe))
    expect_true(all(!is_empty[1:5]))
    expect_true(all(is_empty[-(1:5)]))
    unlink(fn, recursive = TRUE)
})

test_that("Add subset of spots, multiple files in tx spots output", {
    skip_if_not(gdalParquetAvailable())
    fp <- tempdir()
    dir_use <- CosMXOutput(file_path = file.path(fp, "cosmx_test"))

    sfe <- readCosMX(dir_use)
    fn_tx <- formatTxTech(dir_use, tech = "CosMX", z = "all", z_option = "split",
                          file_out = file.path(dir_use, "tx_spots"), return = FALSE)
    sfe <- addSelectTx(sfe, fn_tx, rownames(sfe)[1:5], z_option = "split")
    expect_equal(rowGeometryNames(sfe), c("tx_spots_z0", "tx_spots_z1"))
    z0 <- rowGeometry(sfe, 1L)
    z1 <- rowGeometry(sfe, 2L)
    expect_equal(as.character(st_geometry_type(z0, FALSE)), "MULTIPOINT")
    expect_equal(as.character(st_geometry_type(z1, FALSE)), "MULTIPOINT")
    is_empty0 <- st_is_empty(z0)
    is_empty1 <- st_is_empty(z1)
    expect_false(all(is_empty0))
    expect_false(all(is_empty1))
    expect_true(all(is_empty0[-(1:5)]))
    expect_true(all(is_empty1[-(1:5)]))
    unlink(dir_use, recursive = TRUE)
})

# Final cleanup in case failed test messed with cleanup
fp <- tempdir()
unlink(file.path(fp, "cosmx_test"), recursive = TRUE)
unlink(file.path(fp, "vizgen_test"), recursive = TRUE)
unlink(file.path(fp, "xenium_test"), recursive = TRUE)
