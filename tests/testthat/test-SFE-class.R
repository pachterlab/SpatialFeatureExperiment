library(Matrix)
library(sf)
data("visium_row_col")
# First sample
coords1 <- visium_row_col[visium_row_col$col < 6 & visium_row_col$row < 6, ]
coords1$row <- coords1$row * sqrt(3)
cg <- df2sf(coords1[, c("col", "row")], c("col", "row"), spotDiameter = 0.7)

set.seed(29)
col_inds <- sample(seq_len(13), 13)
row_inds <- sample(seq_len(5), 13, replace = TRUE)
values <- sample(seq_len(5), 13, replace = TRUE)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
colnames(mat) <- coords1$barcode
rownames(mat) <- sample(LETTERS, 5)
rownames(cg) <- colnames(mat)

test_that("Not supplying colGeometries, supplying coordinates in colData", {
    sfe <- SpatialFeatureExperiment(list(counts = mat),
        colData = coords1,
        spatialCoordsNames = c("col", "row")
    )
    expect_equal(colGeometryNames(sfe), "centroids")
    expect_equal(as.character(st_geometry_type(colGeometry(sfe),
        by_geometry = FALSE
    )), "POINT")
})

test_that("Not supplying colGeometries, supplying coordinates in spatialCoords", {
    sfe <- SpatialFeatureExperiment(list(counts = mat),
        spatialCoords = as.matrix(coords1[, c(
            "col",
            "row"
        )])
    )
    expect_equal(colGeometryNames(sfe), "centroids")
    expect_equal(as.character(st_geometry_type(colGeometry(sfe),
        by_geometry = FALSE
    )), "POINT")
})

test_that("Not supplying colGeometries, use spotDiameter for Visium", {
    sfe <- SpatialFeatureExperiment(list(counts = mat),
        colData = coords1,
        spatialCoordsNames = c("col", "row"),
        spotDiameter = 0.7
    )
    expect_equal(colGeometryNames(sfe), "spotPoly")
    expect_equal(as.character(st_geometry_type(colGeometry(sfe),
        by_geometry = FALSE
    )), "POLYGON")
    expect_equal(colGeometry(sfe), cg)
})

test_that("Supplying colGeometries, check centroids", {
    sfe <- SpatialFeatureExperiment(list(counts = mat),
        colGeometries = list(foo = cg)
    )
    expect_equal(colGeometryNames(sfe), "foo")
    centroids <- .sc2cg(spatialCoords(sfe))
    expect_equal(centroids, st_centroid(cg), ignore_attr = TRUE)
})
library(SingleCellExperiment)
test_that("Added object version", {
    sfe <- SpatialFeatureExperiment(list(counts = mat),
                                    colGeometries = list(foo = cg))
    expect_equal(int_metadata(sfe)$SFE_version,
                 packageVersion("SpatialFeatureExperiment"))
})
