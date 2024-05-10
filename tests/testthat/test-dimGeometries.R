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

test_that("rowGeometries setter, one sample", {
    rowGeometries(sfe) <- list(foo = rg_toy)
    rg <- int_elementMetadata(sfe)$rowGeometries$foo
    expect_equal(rg, rg_toy)
})

test_that("rowGeometries setter, one sample, partial replace", {
    foo <- rg_toy[c("d", "j"),]
    rownames(foo) <- c("d", "j")
    rowGeometries(sfe) <- list(foo = foo)
    expect_equal(st_equals(rowGeometry(sfe)[c("d","j"),], foo),
                 list(1, 2), ignore_attr = TRUE)
    expect_true(all(st_is_empty(rowGeometry(sfe)[c(1,3,5),])))

    # The second time, adding more
    bar <- rg_toy[c("v", "o"),]
    rownames(bar) <- c("v", "o")
    rowGeometries(sfe, partial = TRUE) <- list(foo = bar)
    expect_equal(st_equals(rowGeometry(sfe)[2:5,], rg_toy[2:5,]),
                 as.list(1:4), ignore_attr = TRUE)
    expect_true(st_is_empty(rowGeometry(sfe)[1,]))
})

test_that("rowGeometries getter, one sample", {
    internals <- int_elementMetadata(sfe)
    internals[['rowGeometries']] <- make_zero_col_DFrame(nrow(sfe))
    internals[["rowGeometries"]][["foo"]] <- rg_toy
    int_elementMetadata(sfe) <- internals
    rg <- rowGeometries(sfe)
    expect_s4_class(rg, "List")
    expect_equal(names(rg), "foo")
    expect_equal(rg[["foo"]], rg_toy)
})

test_that("rowGeometry setter, one sample", {
    rowGeometry(sfe, "foo") <- rg_toy
    rg <- int_elementMetadata(sfe)$rowGeometries$foo
    expect_equal(rg, rg_toy)
})

test_that("rowGeometry setter, partial replace", {
    foo <- rg_toy[c("d", "j"),]
    rownames(foo) <- c("d", "j")
    rowGeometry(sfe, type = "foo", partial = TRUE) <- foo
    expect_equal(st_equals(rowGeometry(sfe)[c("d","j"),], foo),
                 list(1, 2), ignore_attr = TRUE)
    expect_true(all(st_is_empty(rowGeometry(sfe)[c(1,3,5),])))

    # The second time, adding more
    bar <- rg_toy[c("v", "o"),]
    rownames(bar) <- c("v", "o")
    rowGeometry(sfe, type = "foo", partial = TRUE) <- bar
    expect_equal(st_equals(rowGeometry(sfe)[2:5,], rg_toy[2:5,]),
                 as.list(1:4), ignore_attr = TRUE)
    expect_true(st_is_empty(rowGeometry(sfe)[1,]))
})

test_that("rowGeometry getter, one sample", {
    internals <- int_elementMetadata(sfe)
    internals[['rowGeometries']] <- make_zero_col_DFrame(nrow(sfe))
    internals[["rowGeometries"]][["foo"]] <- rg_toy
    int_elementMetadata(sfe) <- internals
    rg <- rowGeometry(sfe, "foo")
    expect_equal(rg, rg_toy)
})

sfe_visium <- system.file(file.path("extdata", "sfe_visium.rds"),
                          package = "SpatialFeatureExperiment") |>
    readRDS()
rg_toy1 <- cg_toy[1:2,]
rownames(rg_toy1) <- rownames(sfe_visium)
rg_toy2 <- cg_toy[3:4,]
rownames(rg_toy2) <- rownames(sfe_visium)
rg_toy3 <- cg_toy[4:5,]
rownames(rg_toy3) <- rownames(sfe_visium)

test_that("rowGeometries setter, all of multiple samples", {
    rgs <- list(foo = rg_toy1, bar = rg_toy2)
    rowGeometries(sfe_visium, sample_id = "all") <- rgs
    expect_equal(int_elementMetadata(sfe_visium)$rowGeometries |> as.list(),
                 rgs)
})

test_that("rowGeometries setter, some not all of multiple samples", {
    # Subset of samples, first time
    rgs1 <- list(foo = rg_toy1)
    rowGeometries(sfe_visium, sample_id = "sample01") <- rgs1
    expect_equal(rowGeometryNames(sfe_visium), "foo_sample01")
    # Subset of samples, with existing rowGeometries
    rgs2 <- list(bar = rg_toy2)
    rowGeometries(sfe_visium, sample_id = "sample02") <- rgs2
    expect_equal(rowGeometryNames(sfe_visium), c("foo_sample01", "bar_sample02"))
    expect_equal(int_elementMetadata(sfe_visium)$rowGeometries$foo_sample01, rg_toy1)
    expect_equal(int_elementMetadata(sfe_visium)$rowGeometries$bar_sample02, rg_toy2)
})

test_that("rowGeometries setter, partial replace, multiple samples", {
    bar <- rg_toy1[1,,drop = FALSE]
    rownames(bar) <- rownames(sfe_visium)[1]
    baz <- rg_toy2[2,,drop = FALSE]
    rownames(baz) <- rownames(sfe_visium)[2]
    rowGeometries(sfe_visium) <- list(foo_sample01 = bar, foo_sample02 = baz)
    # Second time, adding more
    rowGeometries(sfe_visium, sample_id = "sample01",
                  partial = TRUE) <- list(foo = baz)
    expect_equal(rowGeometry(sfe_visium), rbind(bar, baz), ignore_attr = TRUE)
    rowGeometries(sfe_visium, sample_id = "sample02",
                  partial = TRUE) <- list(foo = bar)
    expect_equal(rowGeometry(sfe_visium, sample_id = "sample02"),
                 rbind(bar, baz), ignore_attr = TRUE)
})

test_that("Set rowGeometries to NULL", {
    rgs <- list(foo = rg_toy1, bar = rg_toy2)
    rowGeometries(sfe_visium, sample_id = "all") <- rgs
    rowGeometries(sfe_visium, sample_id = "all") <- NULL
    expect_equal(length(int_elementMetadata(sfe_visium)$rowGeometries), 0L)
})

test_that("Set rowGeometries of not all samples to NULL", {
    rgs <- list(foo = rg_toy1,
                bar_sample01 = rg_toy2,
                baz_sample02 = rg_toy3)
    rowGeometries(sfe_visium, sample_id = "all") <- rgs
    rowGeometries(sfe_visium, sample_id = "sample02") <- NULL
    expect_equal(rowGeometryNames(sfe_visium), names(rgs)[1:2])
})

test_that("rowGeometries getter, multiple samples", {
    rgs <- list(foo = rg_toy1,
                bar_sample01 = rg_toy2,
                baz_sample02 = rg_toy3)
    rowGeometries(sfe_visium, sample_id = "all") <- rgs
    # All samples
    rgs2 <- rowGeometries(sfe_visium, sample_id = "all")
    expect_equal(rgs2 |> as.list(), rgs)
    # One sample
    rgs3 <- rowGeometries(sfe_visium, sample_id = "sample01")
    expect_equal(as.list(rgs3), rgs[2])
    # Two samples
    rgs4 <- rowGeometries(sfe_visium, sample_id = c("sample01", "sample02"))
    expect_equal(as.list(rgs4), rgs[2:3])
})

test_that("rowGeometry getter, multiple samples", {
    rgs <- list(bar = rg_toy1,
                bar_sample01 = rg_toy2,
                baz_sample02 = rg_toy3)
    rowGeometries(sfe_visium, sample_id = "all") <- rgs
    # Use sample_id = "all" for non-sample-specific rowGeometries
    rg1 <- rowGeometry(sfe_visium, type = "bar", sample_id = "all")
    expect_equal(rg1, rg_toy1)
    # Numeric type
    rg2 <- rowGeometry(sfe_visium, type = 1L, sample_id = "all")
    expect_equal(rg2, rg_toy1)
    expect_error(rowGeometry(sfe_visium, type = 2L, sample_id = "all"),
                 "should not include any sample ID")
    expect_error(rowGeometry(sfe_visium, type = "bar_sample01", sample_id = "all"),
                 "should not include any sample ID")
    # sample-specific
    rg3 <- rowGeometry(sfe_visium, type = "bar", sample_id = "sample01")
    expect_equal(rg3, rg_toy2)
    rg4 <- rowGeometry(sfe_visium, type = 1L, sample_id = "sample01")
    expect_equal(rg4, rg_toy2)
    expect_error(rowGeometry(sfe_visium, type = "baz_sample02", sample_id = "sample01"),
                 "Type does not match sample_id")
})

test_that("rowGeometry setter, multiple samples", {
    rowGeometries(sfe_visium) <- NULL
    # first time
    expect_error(rowGeometry(sfe_visium, type = 1L, sample_id = "sample01") <- rg_toy1,
                 "subscript out of bound for numeric type")
    rowGeometry(sfe_visium, type = "foo", sample_id = "sample01") <- rg_toy1
    expect_equal(rowGeometryNames(sfe_visium), "foo_sample01")
    # later
    rowGeometry(sfe_visium, type = "bar", sample_id = "sample02") <- rg_toy2
    expect_equal(rowGeometryNames(sfe_visium), c("foo_sample01", "bar_sample02"))

    rowGeometries(sfe_visium) <- NULL
    rowGeometry(sfe_visium, type = "foo_sample01", sample_id = "sample01") <- rg_toy1
    expect_equal(rowGeometryNames(sfe_visium), "foo_sample01")
    expect_error(rowGeometry(sfe_visium, type = "foo_sample01", sample_id = "sample02") <- rg_toy1,
                 "Type does not match sample_id")
    expect_error(rowGeometry(sfe_visium, "foo_sample01", sample_id = "all") <- rg_toy1,
                 "should not include any sample ID")
})

test_that("txSpots setter", {
    rowGeometries(sfe) <- NULL
    txSpots(sfe) <- rg_toy
    expect_equal(rowGeometryNames(sfe), "txSpots")
    txSpots(sfe_visium, sample_id = "sample01") <- rg_toy1
    expect_true("txSpots_sample01" %in% rowGeometryNames(sfe_visium))
})

test_that("txSpots getter", {
    txSpots(sfe) <- rg_toy
    expect_equal(txSpots(sfe), rg_toy)
    txSpots(sfe_visium, "sample01") <- rg_toy1
    expect_equal(txSpots(sfe_visium, "sample01"), rg_toy1)
})
