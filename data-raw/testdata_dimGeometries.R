# Toy datasets to unit test dimGeometries getters and setters
library(SpatialExperiment)
library(Matrix)
library(sf)
devtools::load_all()
# Toy SFE object without the geometries. A SPE object without geometries would
# be a valid SFE object.
set.seed(29)
row_inds <- sample(1:5, 5)
col_inds <- sample(1:5, 5)
values <- sample(1:10, 5, replace = TRUE)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
rownames(mat) <- sample(letters, 5)
colnames(mat) <- sample(LETTERS, 5)
coords <- matrix(runif(10), ncol = 2)
spe <- SpatialExperiment(assays = list(counts = mat),
                         sample_id = "sample01",
                         spatialCoords = coords)
sfe <- new("SpatialFeatureExperiment", spe)
saveRDS(sfe, "inst/testdata/sfe_toy.rds")
# Toy colGeometries
cg <- st_sfc(apply(coords, 1, st_point, simplify = FALSE))
cg_sf <- st_sf(geometry = cg, sf_column_name = "geometry")
rownames(cg_sf) <- colnames(mat)
cg_sf2 <- st_buffer(cg_sf, 0.05)
saveRDS(cg_sf, "inst/testdata/cg_toy.rds")
saveRDS(cg_sf2, "inst/testdata/cg_toy2.rds")

# More than one sample_id
spe1 <- SpatialExperiment(assays = list(counts = mat[,1:3]),
                          sample_id = "sample01", spatialCoords = coords[1:3,])
spe2 <- SpatialExperiment(assays = list(counts = mat[,4:5]),
                          sample_id = "sample02", spatialCoords = coords[4:5,])
spe_samples <- cbind(spe1, spe2)
sfe2 <- new("SpatialFeatureExperiment", spe_samples)
saveRDS(sfe2, "inst/testdata/sfe_multi_sample.rds")
