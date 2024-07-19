library(SFEData)
library(sf)

fp <- tempfile()
fn <- XeniumOutput("v2", file_path = fp)

try(sfe <- readXenium(fn))
sfe <- readXenium(fn)
pieces <- readRDS(system.file("extdata/pieces.rds", package = "SpatialFeatureExperiment"))

test_that("Split SFE object by sfc, for one sample", {
    sfes <- splitByCol(sfe, pieces)
    expect_type(sfes, "list")
    expect_length(sfes, 2)
    classes <- vapply(sfes, class, FUN.VALUE = character(1))
    expect_true(all(classes == "SpatialFeatureExperiment"))
    expect_true(all(st_covers(pieces[[1]], colGeometry(sfes[[1]]), sparse = FALSE)))
    expect_true(all(st_covers(pieces[[2]], colGeometry(sfes[[2]]), sparse = FALSE)))
})

# Make 2 samples
sfes <- splitByCol(sfe, pieces)
sfes[[2]] <- changeSampleIDs(sfes[[2]], c(sample01 = "sample02"))
sfe2 <- cbind(sfes[[1]], sfes[[2]])
pieces2 <- readRDS(system.file("extdata/subpieces.rds", package = "SpatialFeatureExperiment"))
pieces2_list <- split(st_geometry(pieces2), pieces2$sample_id)

test_that("Split by list of sfcs, each element for one sample", {
    sfes2 <- splitByCol(sfe2, pieces2_list)
    # there are 4 pieces
    expect_type(sfes2, "list")
    expect_length(sfes2, 4)
    classes <- vapply(sfes2, class, FUN.VALUE = character(1))
    expect_true(all(classes == "SpatialFeatureExperiment"))
    expect_true(all(st_covers(pieces2_list[[1]][[1]], colGeometry(sfes2[[1]]), sparse = FALSE)))
    expect_true(all(st_covers(pieces2_list[[1]][[2]], colGeometry(sfes2[[2]]), sparse = FALSE)))
    expect_true(all(st_covers(pieces2_list[[2]][[1]], colGeometry(sfes2[[3]]), sparse = FALSE)))
    expect_true(all(st_covers(pieces2_list[[2]][[2]], colGeometry(sfes2[[4]]), sparse = FALSE)))
})

test_that("Split by sf, for multiple samples", {
    # The error message
    pieces3 <- pieces2
    names(pieces3)[1] <- "foo"
    expect_error(sfes3 <- splitByCol(sfe2, pieces3), "f must have a column sample_id")
    sfes3 <- splitByCol(sfe2, pieces2)
    expect_type(sfes3, "list")
    expect_length(sfes3, 4)
    classes <- vapply(sfes3, class, FUN.VALUE = character(1))
    expect_true(all(classes == "SpatialFeatureExperiment"))
    expect_true(all(st_covers(st_geometry(pieces2)[[1]], colGeometry(sfes3[[1]]), sparse = FALSE)))
    expect_true(all(st_covers(st_geometry(pieces2)[[4]], colGeometry(sfes3[[2]]), sparse = FALSE)))
    expect_true(all(st_covers(st_geometry(pieces2)[[2]], colGeometry(sfes3[[3]]), sparse = FALSE)))
    expect_true(all(st_covers(st_geometry(pieces2)[[3]], colGeometry(sfes3[[4]]), sparse = FALSE)))
})

test_that("Split different samples into separate SFE objects", {
    sfes4 <- splitSamples(sfe2)
    expect_length(sfes4, 2)
    classes <- vapply(sfes4, class, FUN.VALUE = character(1))
    expect_true(all(classes == "SpatialFeatureExperiment"))
    expect_true(all(st_covers(pieces[[1]], colGeometry(sfes4[[1]]), sparse = FALSE)))
    expect_true(all(st_covers(pieces[[2]], colGeometry(sfes4[[2]]), sparse = FALSE)))
})

cont <- readRDS(system.file("extdata/contiguity.rds", package = "SpatialFeatureExperiment"))
test_that("Split by contiguity of an annotGeometry", {
    cont$sample_id <- "sample01"
    annotGeometry(sfe, "contiguity") <- cont
    sfes5 <- splitContiguity(sfe, annotGeometryName = "contiguity")
    expect_length(sfes5, 2)
    classes <- vapply(sfes5, class, FUN.VALUE = character(1))
    expect_true(all(classes == "SpatialFeatureExperiment"))
    expect_true(all(st_covers(st_geometry(cont)[[1]], colGeometry(sfes5[[1]]), sparse = FALSE)))
    expect_true(all(st_covers(st_union(st_geometry(cont)[2:3]), colGeometry(sfes5[[2]]), sparse = FALSE)))
    # TODO: test when the annotGeometry has elements that are not polygons or multipolygons
    # When it has pieces that are too small
})

unlink(fn, recursive = TRUE)
