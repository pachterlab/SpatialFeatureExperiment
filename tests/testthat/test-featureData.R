library(SFEData)
library(scater)
library(spdep)
library(Voyager)

sfe <- McKellarMuscleData("small")
colGraph(sfe, "visium") <- findVisiumGraph(sfe)

test_that("colFeatureData getter", {
    sfe <- colDataMoransI(sfe, "nCounts")
    cfd <- colFeatureData(sfe)
    expect_s4_class(cfd, "DataFrame")
    expect_equal(names(cfd), c("moran_Vis5A", "K_Vis5A"))
    expect_equal(rownames(cfd), names(colData(sfe)))
    expect_true(is.numeric(cfd["nCounts", "moran_Vis5A"]))
    expect_true(is.na(cfd["nGenes", "moran_Vis5A"]))
})

test_that("rowFeatureData getter", {
    rfd <- rowFeatureData(sfe) # Nothing there yet, don't know what to do with it
    expect_null(rfd)
})

test_that("geometryFeatureData getter", {
    annotGraph(sfe, "myofiber_tri2nb") <-
        findSpatialNeighbors(sfe, type = "myofiber_simplified", MARGIN = 3L,
                             method = "tri2nb", dist_type = "idw",
                             zero.policy = TRUE
        )
    sfe <- annotGeometryMoransI(sfe, features = "area",
                                annotGraphName = "myofiber_tri2nb",
                                annotGeometryName = "myofiber_simplified",
                                zero.policy = TRUE
    )
    gfd <- geometryFeatureData(sfe, "myofiber_simplified", 3L)
    expect_s4_class(gfd, "DataFrame")
    expect_equal(names(gfd), c("moran_Vis5A", "K_Vis5A"))
    expect_equal(rownames(gfd), names(annotGeometry(sfe, "myofiber_simplified")))
    expect_true(is.numeric(gfd["area", "moran_Vis5A"]))
    expect_true(is.na(gfd["perimeter", "moran_Vis5A"]))
})

sfe2 <- sfe[,sfe$in_tissue]
sfe2 <- logNormCounts(sfe2)
sfe2 <- runPCA(sfe2, ncomponents = 10)
annotGraph(sfe2, "myofiber_tri2nb") <-
    findSpatialNeighbors(sfe2, type = "myofiber_simplified", MARGIN = 3L,
                         method = "tri2nb", dist_type = "idw",
                         zero.policy = TRUE
    )
test_that("reducedDimFeatureData getter", {
    sfe2 <- reducedDimMoransI(sfe2, "PCA", components = 1:3)
    rfd <- reducedDimFeatureData(sfe2, "PCA")
    expect_s4_class(rfd, "DataFrame")
    expect_equal(names(rfd), c("moran_Vis5A", "K_Vis5A"))
    expect_equal(rownames(rfd), colnames(reducedDim(sfe2)))
    expect_true(is.numeric(rfd[1:3, "moran_Vis5A"]))
    expect_true(all(is.na(rfd[4:10, "moran_Vis5A"])))
})

test_that("Return NULL when localResults are absent", {
    expect_null(getParams(sfe, "localmoran", local = TRUE,
                          annotGeometryName = "myofiber_simplified"))
    expect_null(getParams(sfe, "localmoran", local = TRUE,
                          colGeometryName = "spotPoly"))
    expect_null(getParams(sfe, "localmoran", local = TRUE))
})

test_that("colFeatureData after deleting columns in colData", {
    sfe <- colDataMoransI(sfe, "nCounts")
    sfe$col <- NULL
    fd <- colFeatureData(sfe)
    expect_equal(rownames(fd), names(colData(sfe)))
})

test_that("colFeatureData after adding columns to colData", {
    sfe <- colDataMoransI(sfe, "nCounts")
    sfe$cluster <- sample(c("foo", "bar", "baz"), ncol(sfe), replace = TRUE)
    fd <- colFeatureData(sfe)
    expect_equal(rownames(fd), names(colData(sfe)))
})

test_that("getParams, gene expression", {
    sfe2 <- runMoransI(sfe2, rownames(sfe2)[1])
    p <- getParams(sfe2, "moran")
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "graph_params"))
})

test_that("getParams, gene expression, local", {
    sfe2 <- runUnivariate(sfe2, "localmoran", rownames(sfe2)[1])
    p <- getParams(sfe2, "localmoran", local = TRUE)
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "p.adjust.method", "graph_params"))
})

test_that("getParams, colData", {
    sfe2 <- colDataMoransI(sfe2, "nCounts")
    p <- getParams(sfe2, "moran", colData = TRUE)
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "graph_params"))
})

test_that("getParams, colData, local", {
    sfe2 <- colDataUnivariate(sfe2, "localmoran", "nCounts")
    p <- getParams(sfe2, "localmoran", local = TRUE, colData = TRUE)
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "p.adjust.method", "graph_params"))
})

test_that("getParams, colGeometry", {
    spotPoly(sfe2)$foo <- rnorm(ncol(sfe2))
    sfe2 <- colGeometryMoransI(sfe2, features = "foo", colGeometryName = "spotPoly")
    p <- getParams(sfe2, "moran", colGeometryName = "spotPoly")
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "graph_params"))
})

test_that("getParams, colGeometry, local", {
    spotPoly(sfe2)$foo <- rnorm(ncol(sfe2))
    sfe2 <- colGeometryUnivariate(sfe2, "localmoran", features = "foo", colGeometryName = "spotPoly")
    p <- getParams(sfe2, "localmoran", local = TRUE, colGeometryName = "spotPoly")
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "p.adjust.method", "graph_params"))
})

test_that("getParams, annotGeometry", {
    sfe2 <- annotGeometryMoransI(sfe2, features = "area",
                                 annotGraphName = "myofiber_tri2nb",
                                 annotGeometryName = "myofiber_simplified",
                                 zero.policy = TRUE
    )
    p <- getParams(sfe2, "moran", annotGeometryName = "myofiber_simplified")
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "graph_params"))
})

test_that("getParams, annotGeometry, local", {
    sfe2 <- annotGeometryUnivariate(sfe2, "localmoran", features = "area",
                                    annotGraphName = "myofiber_tri2nb",
                                    annotGeometryName = "myofiber_simplified",
                                    zero.policy = TRUE
    )
    p <- getParams(sfe2, "localmoran", local = TRUE, annotGeometryName = "myofiber_simplified")
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "p.adjust.method", "graph_params"))
})

test_that("getParams, reducedDim", {
    sfe2 <- reducedDimMoransI(sfe2, "PCA", components = 1:3)
    p <- getParams(sfe2, "moran", reducedDimName = "PCA")
    expect_named(p, c("name", "package", "version", "zero.policy", "include_self",
                      "graph_params"))
})

sfe3 <- McKellarMuscleData("small2")
sfe3 <- sfe3[,sfe3$in_tissue]
colGraph(sfe3, "visium") <- findVisiumGraph(sfe3, style = "S")
sfe3 <- logNormCounts(sfe3)
sfe3 <- runPCA(sfe3, ncomponents = 10)
annotGraph(sfe3, "myofiber_tri2nb") <-
    findSpatialNeighbors(sfe3, type = "myofiber_simplified", MARGIN = 3L,
                         method = "tri2nb", dist_type = "idw",
                         zero.policy = TRUE
    )
test_that("Deal with featureData in cbind, when both SFE objects have featureData", {
    sfe2 <- runMoransI(sfe2, rownames(sfe2)[1])
    sfe3 <- runMoransI(sfe3, rownames(sfe3)[1])
    sfe2 <- colDataMoransI(sfe2, "nCounts")
    sfe3 <- colDataMoransI(sfe3, "nCounts")
    sfe2 <- annotGeometryMoransI(sfe2, features = "area",
                                 annotGraphName = "myofiber_tri2nb",
                                 annotGeometryName = "myofiber_simplified",
                                 zero.policy = TRUE
    )
    sfe3 <- annotGeometryMoransI(sfe3, features = "area",
                                 annotGraphName = "myofiber_tri2nb",
                                 annotGeometryName = "myofiber_simplified",
                                 zero.policy = TRUE
    )
    sfe2 <- reducedDimMoransI(sfe2, "PCA", components = 1:3)
    sfe3 <- reducedDimMoransI(sfe3, "PCA", components = 1:3)
    sfe_c <- cbind(sfe2, sfe3)

    name_expect <- c("moran_Vis5A", "K_Vis5A", "moran_sample02", "K_sample02")
    cfd <- colFeatureData(sfe_c)
    expect_equal(names(cfd), name_expect)
    afd <- geometryFeatureData(sfe_c, "myofiber_simplified", 3L)
    expect_equal(names(afd), name_expect)
    rfd <- reducedDimFeatureData(sfe_c, "PCA")
    expect_equal(names(rfd), name_expect)
})

test_that("Deal with featureData in cbind, when only one SFE objects has featureData", {
    sfe2 <- runMoransI(sfe2, rownames(sfe2)[1])
    sfe2 <- colDataMoransI(sfe2, "nCounts")
    sfe2 <- annotGeometryMoransI(sfe2, features = "area",
                                 annotGraphName = "myofiber_tri2nb",
                                 annotGeometryName = "myofiber_simplified",
                                 zero.policy = TRUE
    )
    sfe2 <- reducedDimMoransI(sfe2, "PCA", components = 1:3)
    sfe_c <- cbind(sfe2, sfe3)
    name_expect <- c("moran_Vis5A", "K_Vis5A")
    cfd <- colFeatureData(sfe_c)
    expect_equal(names(cfd), name_expect)
    afd <- geometryFeatureData(sfe_c, "myofiber_simplified", 3L)
    expect_equal(names(afd), name_expect)
    rfd <- reducedDimFeatureData(sfe_c, "PCA")
    expect_equal(names(rfd), name_expect)
})
