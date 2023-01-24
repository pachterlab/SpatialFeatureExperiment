outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
bc_flou1 <- read.csv(file.path(outdir, "sample01", "outs", "spatial", 
                               "barcode_fluorescence_intensity.csv"))
sp_enr1 <- read.csv(file.path(outdir, "sample01", "outs", "spatial", 
                              "spatial_enrichment.csv"))
pos1 <- read.csv(file.path(outdir, "sample01", "outs", "spatial", 
                           "tissue_positions.csv"))

bc_flou2 <- read.csv(file.path(outdir, "sample02", "outs", "spatial", 
                               "barcode_fluorescence_intensity.csv"))
sp_enr2 <- read.csv(file.path(outdir, "sample02", "outs", "spatial", 
                              "spatial_enrichment.csv"))
pos2 <- read.csv(file.path(outdir, "sample02", "outs", "spatial", 
                           "tissue_positions.csv"))

samples <- file.path(outdir, paste0("sample0", 1:2))

rd1 <- sp_enr1[,4:9]
rd2 <- sp_enr2[,4:9]
names(rd1) <- paste(names(rd1), "sample01", sep = "_")
names(rd2) <- paste(names(rd2), "sample02", sep = "_")

rd_expect <- cbind(Feature.Type = sp_enr1$Feature.Type, rd1, rd2)

cd_expect <- rbind(bc_flou1, bc_flou2)[, 3:8]

test_that("Correctly read Space Ranger output", {
    sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
    # Very basic one
    expect_s4_class(sfe, "SpatialFeatureExperiment")
    expect_equal(sampleIDs(sfe), c("sample01", "sample02"))
    expect_equal(colGeometryNames(sfe), "spotPoly")
    expect_equal(colGraphNames(sfe, "sample01"), "visium")
    expect_equal(colGraphNames(sfe, "sample02"), "visium")
    expect_equal(as.data.frame(rowData(sfe)[,-1]), rd_expect,
                 ignore_attr = "row.names")
    expect_equal(as.data.frame(colData(sfe)[,5:10]), cd_expect,
                 ignore_attr = "row.names")
})
