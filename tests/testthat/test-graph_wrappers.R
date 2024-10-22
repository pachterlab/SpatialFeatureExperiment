library(spdep)
set.SubgraphOption(FALSE) # I don't care in this case
# Unit test spdep and Visium graph wrappers
sfe2 <- readRDS(system.file("extdata/sfe_multi_sample.rds",
    package = "SpatialFeatureExperiment"
))
cgr1 <- readRDS(system.file("extdata/colgraph1.rds",
    package = "SpatialFeatureExperiment"
))

test_that("Get the correct graph and attr for reconstruction", {
    g <- findSpatialNeighbors(sfe2, sample_id = "sample01",
        type = "spatialCoords",
        MARGIN = 2, method = "tri2nb"
    )
    expect_equal(g, cgr1, ignore_attr = TRUE)
    attrs_reconst <- attr(g, "method")
    expect_equal(names(attrs_reconst), c("FUN", "package", "args"))
    expect_equal(attrs_reconst$FUN, "findSpatialNeighbors")
    expect_equal(attrs_reconst$package[[1]], "SpatialFeatureExperiment")
    expect_equal(attrs_reconst$package[[2]], packageVersion("SpatialFeatureExperiment"))
    expect_equal(attrs_reconst$args$dist_type, "none")
    expect_equal(attrs_reconst$args$style, "W")
    expect_true("row.names" %in% names(attrs_reconst$args))
    expect_equal(attrs_reconst$args$method, "tri2nb")
})

test_that("Use distance weighting", {
    g <- findSpatialNeighbors(sfe2, "sample01",
        type = "spatialCoords",
        MARGIN = 2, method = "tri2nb",
        dist_type = "idw"
    )
    expect_equal(g$neighbours, cgr1$neighbours, ignore_attr = TRUE)
    expect_false(isTRUE(all.equal(g$weights, cgr1$weights, check.attributes = FALSE)))
    attrs_reconst <- attr(g, "method")
    expect_equal(attrs_reconst$args$dist_type, "idw")
    expect_equal(attrs_reconst$args$style, "raw")
})

library(SFEData)
sfe_muscle1 <- McKellarMuscleData("small")
sfe_muscle2 <- McKellarMuscleData("small2")
sfe <- cbind(sfe_muscle1, sfe_muscle2)

dist_types <- c("none", "idw", "exp", "dpd")
styles <- c("raw", "W", "B", "C", "U", "minmax", "S")
test_that("Exact Bioc methods for knn return same results as spdep methods", {
    # Types of distance weights
    for (d in dist_types) {
        s <- "W"
        cat("Testing dist_type", d, "style", s, "\n")
        g1 <- findSpatialNeighbors(sfe, k = 3, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "knearneigh",
                                   nn_method = "spdep", dmax = 200)
        g2 <- findSpatialNeighbors(sfe, k = 3, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "knearneigh",
                                   nn_method = "bioc", dmax = 200)
        expect_equal(g1, g2, ignore_attr = TRUE)
    }
    # Types of normalization
    for (s in styles) {
        d <- "idw"
        cat("Testing dist_type", d, "style", s, "\n")
        g1 <- findSpatialNeighbors(sfe, k = 3, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "knearneigh",
                                   nn_method = "spdep")
        g2 <- findSpatialNeighbors(sfe, k = 3, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "knearneigh",
                                   nn_method = "bioc")
        expect_equal(g1, g2, ignore_attr = TRUE)
    }
})

test_that("Exact Bioc methods for dnearneigh return same results as spdep methods", {
    # Types of distance weights
    for (d in dist_types) {
        s <- "W"
        cat("Testing dist_type", d, "style", s, "\n")
        g1 <- findSpatialNeighbors(sfe, d1 = 50, d2 = 150, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "dnearneigh",
                                   nn_method = "spdep", dmax = 200)
        g2 <- findSpatialNeighbors(sfe, d1 = 50, d2 = 150, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "dnearneigh",
                                   nn_method = "bioc", dmax = 200)
        expect_equal(g1, g2, ignore_attr = TRUE)
    }
    # Types of normalization
    for (s in styles) {
        d <- "idw"
        cat("Testing dist_type", d, "style", s, "\n")
        g1 <- findSpatialNeighbors(sfe, d1 = 50, d2 = 150, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "dnearneigh",
                                   nn_method = "spdep")
        g2 <- findSpatialNeighbors(sfe, d1 = 50, d2 = 150, sample_id = "all", MARGIN = 3,
                                   dist_type = d, style = s,
                                   type = "myofiber_simplified",
                                   method = "dnearneigh",
                                   nn_method = "bioc")
        expect_equal(g1, g2, ignore_attr = TRUE)
    }
})

test_that("Error when dmax is not specified for DPD", {
    expect_error(g <- findSpatialNeighbors(sfe, dist_type = "dpd", d2 = 150,
                                           sample_id = "all", MARGIN = 3,
                                           type = "myofiber_simplified",
                                           method = "dnearneigh",
                                           nn_method = "bioc"),
                 "DPD weights require a maximum distance threshold")
    expect_error(g <- findSpatialNeighbors(sfe, dist_type = "dpd", d2 = 150,
                                           sample_id = "all", MARGIN = 3,
                                           type = "myofiber_simplified",
                                           method = "dnearneigh",
                                           nn_method = "bioc", dmax = -2),
                 "DPD weights require a positive")
    expect_error(g <- findSpatialNeighbors(sfe, dist_type = "dpd", k = 10,
                                           sample_id = "all", MARGIN = 3,
                                           type = "myofiber_simplified",
                                           method = "knearneigh",
                                           nn_method = "bioc"),
                 "DPD weights require a maximum distance threshold")
    expect_error(g <- findSpatialNeighbors(sfe, dist_type = "dpd", k = 10,
                                           sample_id = "all", MARGIN = 3,
                                           type = "myofiber_simplified",
                                           method = "knearneigh",
                                           nn_method = "bioc", dmax = -2),
                 "DPD weights require a positive")
})

sfe_visium <- readRDS(system.file("extdata/sfe_visium.rds",
    package = "SpatialFeatureExperiment"
))
g_visium <- readRDS(system.file("extdata/colgraph_visium.rds",
    package = "SpatialFeatureExperiment"
))
test_that("Correct Visium graph", {
    g <- findVisiumGraph(sfe_visium, "sample01")
    expect_equal(g, g_visium, ignore_attr = TRUE)
    attrs_reconst <- attr(g, "method")
    expect_equal(attrs_reconst$FUN, "findVisiumGraph")
    expect_equal(attrs_reconst$args$style, "W")
})

test_that("Correct Visium HD graph", {
    testthat::skip_on_ci()
    dir <- "~/WoundAnalysis/Visium-HD data/YVW01_binned_outputs/"
    sfe <- readVisiumHD(dir, bin_size = 16)
    g <- findVisiumHDGraph(sfe)
    expect_s3_class(g, "listw")
})
