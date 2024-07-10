# Toy Visium SFE dataset to unit test Visium specific functions
library(tidyverse)
library(spdep)
library(Matrix)
library(SpatialExperiment)
devtools::load_all()
data("visium_row_col")
coords <- visium_row_col %>%
  filter(col < 6, row < 6)
coords_mat <- as.matrix(coords[,c("col", "row")])
# make hexagonal grid
coords_mat[,"row"] <- coords_mat[,"row"] * sqrt(3)
# Make two samples, physical rows 1 and 2 go to the first one
sample01_ind <- coords_mat[,"row"] < 2.5 * sqrt(3)
g1 <- nb2listw(dnearneigh(coords_mat[sample01_ind,], 1.9, 2.1,
                          row.names = coords$barcode[sample01_ind]))
g2 <- nb2listw(dnearneigh(coords_mat[!sample01_ind,], 1.9, 2.1))
attr(g1, "method") <- list(FUN = "findVisiumGraph",
                           package = "SpatialFeatureExperiment",
                           args = list(barcode_allow_list = NULL,
                                       style = "W",
                                       zero.policy = NULL,
                                       sample_id = "sample01"))
# g2 doesn't have the attributes on purpose for unit test purposes
saveRDS(g1, "inst/extdata/colgraph_visium.rds")
saveRDS(g2, "inst/extdata/colgraph_visium2.rds")

set.seed(29)
col_inds <- sample(seq_len(13), 5)
row_inds <- sample(seq_len(2), 5, replace = TRUE)
values <- sample(seq_len(10), 5)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
rownames(mat) <- sample(LETTERS, 2)
colnames(mat) <- coords$barcode
spe1 <- SpatialExperiment(assays = list(counts = mat[,sample01_ind]),
                         sample_id = "sample01",
                         spatialCoords = coords_mat[sample01_ind,])
spe2 <- SpatialExperiment(assays = list(counts = mat[,!sample01_ind]),
                          sample_id = "sample02",
                          spatialCoords = coords_mat[!sample01_ind,])
spe <- cbind(spe1, spe2)
sfe <- new("SpatialFeatureExperiment", spe)
saveRDS(sfe, "inst/extdata/sfe_visium.rds")

# Graph when one vertex is removed
g1_sub <- nb2listw(dnearneigh(spatialCoords(spe1)[-1,], 1.9, 2.1,
                              row.names = colnames(spe1)[-1]))
attr(g1_sub, "method") <- list(FUN = "findVisiumGraph",
                               package = "SpatialFeatureExperiment",
                               args = list(barcode_allow_list = NULL,
                                           style = "W",
                                           zero.policy = NULL,
                                           sample_id = "sample01"))
saveRDS(g1_sub, "inst/extdata/colgraph_visium_sub.rds")
