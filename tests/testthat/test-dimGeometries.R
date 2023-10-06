# Unit test dimGeometry/ies getters and setters
library(SingleCellExperiment)
library(S4Vectors)
library(sf)
sfe <- readRDS(system.file("extdata/sfe_toy.rds",
    package = "SpatialFeatureExperiment"
))

test_that("Get List of length 0 when dimGeometries are absent", {
    foo <- dimGeometries(sfe, 2)
    expect_true(is(foo, "List"))
    expect_equal(length(foo), 0L)
})

cg_toy <- readRDS(system.file("extdata/cg_toy.rds",
    package = "SpatialFeatureExperiment"
))
cg_toy2 <- readRDS(system.file("extdata/cg_toy2.rds",
    package = "SpatialFeatureExperiment"
))

test_that("colGeometries setter", {
    colGeometries(sfe) <- list(coords = cg_toy, buffered = cg_toy2)
    cgs <- int_colData(sfe)$colGeometries
    expect_true(is(cgs, "DFrame"))
    expect_equal(names(cgs), c("coords", "buffered"))
    expect_true(is(cgs$coords, "sf"))
    expect_true(is(cgs$buffered, "sf"))
})

test_that("colGeometry setter", {
    colGeometry(sfe, "coords") <- cg_toy
    cgs <- int_colData(sfe)$colGeometries
    expect_true(is(cgs, "DFrame"))
    expect_equal(names(cgs), "coords")
    expect_true(is(cgs$coords, "sf"))
})

sfe2 <- sfe
int_colData(sfe2)$colGeometries <- make_zero_col_DFrame(nrow = ncol(sfe2))
int_colData(sfe2)$colGeometries$coords <- cg_toy
int_colData(sfe2)$colGeometries$buffered <- cg_toy2

test_that("colGeometries getter", {
    cgs <- colGeometries(sfe2)
    expect_true(is(cgs, "List"))
    expect_equal(length(cgs), 2L)
    expect_equal(names(cgs), c("coords", "buffered"))
})

test_that("colGeometry getter", {
    cg <- colGeometry(sfe2)
    expect_equal(cg, cg_toy)
    cg2 <- colGeometry(sfe2, "buffered")
    expect_equal(cg2, cg_toy2)
    cg3 <- colGeometry(sfe2, 2L)
    expect_equal(cg3, cg_toy2)
    expect_error(colGeometry(sfe2, "purr"), "not in")
})

test_that("colGeometryNames", {
    nms <- colGeometryNames(sfe2)
    expect_equal(nms, c("coords", "buffered"))
    colGeometryNames(sfe2) <- c("foo", "bar")
    expect_equal(names(int_colData(sfe2)$colGeometries), c("foo", "bar"))
})

# More than one sample_id
sfe3 <- readRDS(system.file("extdata/sfe_multi_sample.rds",
    package = "SpatialFeatureExperiment"
))

test_that("colGeometry setter for one of the two samples (not already present)", {
    # colGeometry not already present
    colGeometry(sfe3,
        type = "coords", sample_id = "sample01",
        withDimnames = FALSE
    ) <- cg_toy[seq_len(3), ]
    expect_equal(
        int_colData(sfe3)$colGeometries$coords[seq_len(3), ],
        cg_toy[seq_len(3), ]
    )
    expect_true(all(st_is_empty(int_colData(sfe3)$colGeometries$coords[4:5, ])))
})

test_that("colGeometry setter for one of the two samples (already present)", {
    # colGeometry already present
    sfe3 <- addVisiumSpotPoly(sfe3, 0.3)
    cg_orig <- int_colData(sfe3)$colGeometries$spotPoly
    colGeometry(sfe3, "spotPoly", sample_id = "sample01") <- cg_toy[seq_len(3), ]
    expect_equal(int_colData(sfe3)$colGeometries$spotPoly[seq_len(3), ],
        cg_toy[seq_len(3), ])
    expect_equal(int_colData(sfe3)$colGeometries$spotPoly[4:5, ], cg_orig[4:5, ])
})

test_that("colGeometry getter for one of the two samples", {
    int_colData(sfe3)$colGeometries <- make_zero_col_DFrame(nrow = ncol(sfe3))
    int_colData(sfe3)$colGeometries$coords <- cg_toy
    coords_sample02 <- colGeometry(sfe3, "coords", "sample02")
    expect_equal(nrow(coords_sample02), 2L)
    cg_expect <- cg_toy[4:5, ]
    rownames(cg_expect) <- rownames(cg_toy)[4:5]
    expect_equal(coords_sample02, cg_expect)
})

test_that("When rownames of value don't match those of sfe", {
    cg_sub <- cg_mut <- cg_toy[1:3,]
    rownames(cg_sub) <- rownames(cg_toy)[1:3]
    rownames(cg_mut) <- c(rownames(cg_toy)[1:2], "D")
    # It's not a subset, get correct error message
    expect_error(colGeometry(sfe, "foo") <- cg_mut,
                 "should all be in")
    # It's the same set, in different orders, do the reordering
    rns <- sample(rownames(cg_toy), 5, replace = FALSE)
    cg_perm <- cg_toy[rns,]
    rownames(cg_perm) <- rns
    colGeometry(sfe, "foo") <- cg_perm
    expect_equal(colGeometry(sfe, "foo"), cg_toy)
    # It's a subset, do the reordering
    cg_sub_perm <- cg_sub[c(2,3,1),]
    rownames(cg_sub_perm) <- rownames(cg_sub)[c(2,3,1)]
    colGeometry(sfe, "bar") <- cg_sub_perm
    expect_equal(colGeometry(sfe, "bar")[1:3,], cg_sub, ignore_attr = "row.names")
})

rg_toy <- cg_toy
rownames(rg_toy) <- rownames(sfe)

test_that("rowGeometries setter", {
    rowGeometries(sfe) <- list(foo = rg_toy)
    rg <- int_elementMetadata(sfe)$rowGeometries$foo
    expect_equal(rg, rg_toy)
})

test_that("rowGeometries getter", {
    internals <- int_elementMetadata(sfe)
    internals[['rowGeometries']] <- make_zero_col_DFrame(nrow(sfe))
    internals[["rowGeometries"]][["foo"]] <- rg_toy
    int_elementMetadata(sfe) <- internals
    rg <- rowGeometries(sfe)
    expect_s4_class(rg, "List")
    expect_equal(names(rg), "foo")
    expect_equal(rg[["foo"]], rg_toy)
})

test_that("rowGeometry setter", {
    rowGeometry(sfe, "foo") <- rg_toy
    rg <- int_elementMetadata(sfe)$rowGeometries$foo
    expect_equal(rg, rg_toy)
})

test_that("rowGeometry getter", {
    internals <- int_elementMetadata(sfe)
    internals[['rowGeometries']] <- make_zero_col_DFrame(nrow(sfe))
    internals[["rowGeometries"]][["foo"]] <- rg_toy
    int_elementMetadata(sfe) <- internals
    rg <- rowGeometry(sfe, "foo")
    expect_equal(rg, rg_toy)
})
