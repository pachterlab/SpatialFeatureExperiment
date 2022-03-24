# Toy Visium SFE dataset to unit test Visium specific functions
library(tidyverse)
library(spdep)
library(Matrix)
devtools::load_all()
data("visium_row_col")
coords <- visium_row_col %>%
  filter(col < 6, row < 6)
coords_mat <- as.matrix(coords[,c("col", "row")])
# make hexagonal grid
coords_mat[,"row"] <- coords_mat[,"row"] * sqrt(3)
g <- dnearneigh(coords_mat, 1.9, 2.1)
g_listw <- nb2listw(g)
saveRDS(g_listw, "inst/testdata/colgraph_visium.rds")

set.seed(29)
col_inds <- sample(1:13, 5)
row_inds <- sample(1:2, 5, replace = TRUE)
values <- sample(1:10, 5)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
rownames(mat) <- sample(LETTERS, 2)
colnames(mat) <- coords$barcode
spe <- SpatialExperiment(assays = list(counts = mat),
                         sample_id = "sample01",
                         spatialCoords = coords_mat)
sfe <- new("SpatialFeatureExperiment", spe)
saveRDS(sfe, "inst/testdata/sfe_visium.rds")
