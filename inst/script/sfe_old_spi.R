# Make rds file of sfe object with images but in old SpatRasterImage
# Use 3.18 release version here
devtools::load_all()
fp <- system.file("extdata/sample01", package = "SpatialFeatureExperiment")
sfe <- read10xVisiumSFE(fp, type = "sparse")
saveRDS(sfe, "inst/extdata/sfe_old_spi.rds")
