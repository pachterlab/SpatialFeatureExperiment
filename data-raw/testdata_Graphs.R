# Toy datasets to unit test spatialGraph(s) functions
library(spdep)
library(SpatialExperiment)
devtools::load_all()
# See the data-raw/testdata_Geometries.R file
sfe2 <- readRDS("inst/testdata/sfe_multi_sample.rds")
coords <- spatialCoords(sfe2)
# Two samples, colGraphs
cgr1 <- nb2listw(tri2nb(coords[1:3,], row.names = colnames(sfe2)[1:3]))
cgr2 <- nb2listw(knn2nb(knearneigh(coords[4:5,]), row.names = colnames(sfe2)[4:5]))
saveRDS(cgr1, "inst/testdata/colgraph1.rds")
saveRDS(cgr2, "inst/testdata/colgraph2.rds")
# Some random annotGraphs
set.seed(89)
coords <- matrix(runif(12), ncol = 2)
agr1 <- nb2listw(tri2nb(coords[1:3,]))
agr2 <- nb2listw(tri2nb(coords[4:6,]))
saveRDS(agr1, "inst/testdata/annotgraph1.rds")
saveRDS(agr2, "inst/testdata/annotgraph2.rds")
