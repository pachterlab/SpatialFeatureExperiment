# Toy datasets to unit test *Geometries getters and setters
library(SpatialExperiment)
library(Matrix)
library(sf)
devtools::load_all()
# Toy SFE object without the geometries. A SPE object without geometries would
# be a valid SFE object.
set.seed(29)
row_inds <- sample(seq_len(5), 5)
col_inds <- sample(seq_len(5), 5)
values <- sample(seq_len(10), 5, replace = TRUE)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
rownames(mat) <- sample(letters, 5)
colnames(mat) <- sample(LETTERS, 5)
coords <- matrix(runif(10), ncol = 2)
spe <- SpatialExperiment(assays = list(counts = mat),
                         sample_id = "sample01",
                         spatialCoords = coords)
sfe <- new("SpatialFeatureExperiment", spe)
saveRDS(sfe, "inst/extdata/sfe_toy.rds")
# Toy colGeometries
cg <- st_sfc(apply(coords, 1, st_point, simplify = FALSE))
cg_sf <- st_sf(geometry = cg, sf_column_name = "geometry")
rownames(cg_sf) <- colnames(mat)
cg_sf2 <- st_buffer(cg_sf, 0.05)
saveRDS(cg_sf, "inst/extdata/cg_toy.rds")
saveRDS(cg_sf2, "inst/extdata/cg_toy2.rds")

# More than one sample_id
spe1 <- SpatialExperiment(assays = list(counts = mat[,seq_len(3)]),
                          sample_id = "sample01",
                          spatialCoords = coords[seq_len(3),])
spe2 <- SpatialExperiment(assays = list(counts = mat[,4:5]),
                          sample_id = "sample02", spatialCoords = coords[4:5,])
spe_samples <- cbind(spe1, spe2)
sfe2 <- new("SpatialFeatureExperiment", spe_samples)
saveRDS(sfe2, "inst/extdata/sfe_multi_sample.rds")

# annotGeometries
# The difference is that the number of rows is not regulated
ag <- st_sf(geometry = st_convex_hull(st_combine(cg_sf)),
            sample_id = "sample01",
            sf_column_name = "geometry")
# For different sample_ids
ag1 <- st_cast(st_combine(cg_sf[seq_len(3),]), "LINESTRING")
ag2 <- st_cast(st_combine(cg_sf[4:5,]), "LINESTRING")
ag_samples <- st_sf(geometry = c(ag1, ag2),
                    sample_id = c("sample01", "sample02"),
                    sf_column_name = "geometry")
saveRDS(ag, "inst/extdata/ag.rds")
saveRDS(ag_samples, "inst/extdata/ag_samples.rds")
