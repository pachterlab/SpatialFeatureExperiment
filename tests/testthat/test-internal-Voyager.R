library(SFEData)

sfe <- McKellarMuscleData("small")
# To make sure there's no code lingering around stuck to old "symbol" requirement
names(rowData(sfe))[names(rowData(sfe)) == "symbol"] <- "Symbol"

spotPoly(sfe)$foo <- sfe$nCounts
gene_use <- rownames(sfe)[1]

test_that(".check_features gives correct results", {
    res <- .check_features(sfe, c("nCounts", gene_use, "foo"),
                           colGeometryName = "spotPoly")
    expect_type(res, "list")
    expect_named(res, c("assay", "coldata", "colgeom"))
    expect_equal(res$assay, gene_use)
    expect_equal(res$coldata, "nCounts")
    expect_equal(res$colgeom, "foo")
    expect_error(.check_features(sfe, "bar"),
                 "None of the features")
})

symbol1 <- rowData(sfe)$Symbol[1]
dup_ind <- which(duplicated(rowData(sfe)$Symbol))
symbol2 <- rowData(sfe)$Symbol[dup_ind[1]]

test_that(".check_features use swap_rownames and gene symbol input", {
    res <- .check_features(sfe, symbol1, swap_rownames = "Symbol")
    expect_equal(res$assay, gene_use)
    expect_true(all(vapply(res[c("coldata", "colgeom")], function(x) length(x) == 0L,
                           FUN.VALUE = logical(1))))
    # When there're duplicate symbols
    expect_warning(res <- .check_features(sfe, symbol2, swap_rownames = "Symbol"),
                   "duplicated")
    expect_equal(res$assay, rownames(sfe)[rowData(sfe)$Symbol == symbol2][1])
})

test_that(".symbol2id gives correct results", {
    res <- .symbol2id(sfe, symbol1, swap_rownames = "Symbol")
    expect_equal(res, gene_use)
    expect_warning(res <- .symbol2id(sfe, symbol2, "Symbol"),
                   "duplicated")
    expect_equal(res, rownames(sfe)[rowData(sfe)$Symbol == symbol2][1])
})
