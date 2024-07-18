sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
    package = "SpatialFeatureExperiment"
))
g_visium <- readRDS(system.file("extdata/colgraph_visium.rds",
    package = "SpatialFeatureExperiment"
))
no_valid_barcodes <- "Invalid barcodes removed. Valid barcodes: 0/5"
one_valid_barcodes <- "Invalid barcodes removed. Valid barcodes: 1/5"

check_sfe_graph <- function(actual_g, expected_g) {
    expect_equal(actual_g, expected_g, ignore_attr = TRUE)
    attrs_reconst <- attr(actual_g, "method")
    expect_equal(attrs_reconst$FUN, "findVisiumGraph")
    expect_equal(attrs_reconst$args$style, "W")
}

test_that("Correct Visium graph with all valid barcodes", {
    g <- findVisiumGraph(sfe_visium, "sample01")
    check_sfe_graph(g, g_visium)
})

test_that("Visium graph with all incorrect barcodes", {
    data(visium_row_col_v5)
    error <- "After filtering by valid barcode, there were none left, try passing a different barcode_allow_list."
    expect_error(
        expect_warning(
            findVisiumGraph(sfe_visium, "sample01", visium_row_col_v5), no_valid_barcodes
        ), error
    )
})

test_that("Visium graph with one correct barcode", {
    data(visium_row_col_v5)
    changed_sfe <- sfe_visium
    colnames(changed_sfe)[[1]] <- "AACAATCCGAGTGGAC-1"
    # Expects error from dnearneigh but as it's a 3rd part library just checking for an error - not message.
    expect_error(
        expect_warning(
            findVisiumGraph(changed_sfe, "sample01", visium_row_col_v5), one_valid_barcodes
        )
    )
})
