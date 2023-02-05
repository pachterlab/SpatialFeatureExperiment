# Make my own toy 10x dataset
# With spatial enrichment file
library(rjson)
library(Matrix)
library(tidyverse)
library(EBImage)

if (!file.exists("visium_prostate.tar.gz"))
  download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_filtered_feature_bc_matrix.tar.gz", 
                destfile = "visium_prostate.tar.gz")
if (!file.exists("visium_prostate_spatial.tar.gz"))
  download.file("https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/Visium_FFPE_Human_Prostate_IF/Visium_FFPE_Human_Prostate_IF_spatial.tar.gz", 
                destfile = "visium_prostate_spatial.tar.gz")
if (!dir.exists("outs")) {
  dir.create("outs")
  system("tar -xvf visium_prostate.tar.gz -C outs")
  system("tar -xvf visium_prostate_spatial.tar.gz -C outs")
}

# Subset colData and rowData
bc_fluo <- read.csv("outs/spatial/barcode_fluorescence_intensity.csv")
sp_enr <- read.csv("outs/spatial/spatial_enrichment.csv")
tissue_pos <- read.csv("outs/spatial/tissue_positions.csv")

tissue_pos2 <- tissue_pos |> 
  filter(between(array_col, 101, 105),
         between(array_row, 11, 15), # sample02 is 11 to 15
         in_tissue > 0)
outpath <- "inst/extdata/sample02/outs"
bc_fluo2 <- bc_fluo[match(tissue_pos2$barcode, bc_fluo$barcode),]
sp_enr2 <- sp_enr[1:5,]

dir.create(outpath)
dir.create(file.path(outpath, "spatial"))

write.csv(bc_fluo2, file.path(outpath, "spatial", "barcode_fluorescence_intensity.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(sp_enr2, file.path(outpath, "spatial", "spatial_enrichment.csv"),
          quote = FALSE, row.names = FALSE)
write.csv(tissue_pos2, file.path(outpath, "spatial", "tissue_positions.csv"),
          quote = FALSE, row.names = FALSE)

# Subset matrix
bcs <- read.table("outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
features <- read.table("outs/filtered_feature_bc_matrix/features.tsv.gz", sep = "\t")
mat <- readMM("outs/filtered_feature_bc_matrix/matrix.mtx.gz")

colnames(mat) <- bcs$V1
rownames(mat) <- features$V1

features2 <- features[match(sp_enr2$Feature.ID, features$V1),]
mat2 <- mat[features2$V1, bc_fluo2$barcode]

dir.create(file.path(outpath, "filtered_feature_bc_matrix"))
writeLines(colnames(mat2), file(file.path(outpath, "filtered_feature_bc_matrix", "barcodes.tsv")))
write.table(features2, file.path(outpath, "filtered_feature_bc_matrix", "features.tsv"),
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
writeMM(mat2, file.path(outpath, "filtered_feature_bc_matrix", "matrix.mtx"))

# Crop images
scalefactors <- fromJSON(file = "outs/spatial/scalefactors_json.json")
fullres_bbox <- c(range(tissue_pos2$pxl_col_in_fullres), range(tissue_pos2$pxl_row_in_fullres))
names(fullres_bbox) <- c("xmin", "xmax", "ymin", "ymax")
lowres_bbox <- round(fullres_bbox * scalefactors$tissue_lowres_scalef)
hires_bbox <- round(fullres_bbox * scalefactors$tissue_hires_scalef)

lowres <- readImage("outs/spatial/tissue_lowres_image.png")
hires <- readImage("outs/spatial/tissue_hires_image.png")
aligned <- readImage("outs/spatial/aligned_fiducials.jpg")
detected <- readImage("outs/spatial/detected_tissue_image.jpg")

lowres2 <- lowres[lowres_bbox["xmin"]:lowres_bbox["xmax"], 
                  lowres_bbox["ymin"]:lowres_bbox["ymax"],]
hires2 <- hires[hires_bbox["xmin"]:hires_bbox["xmax"], 
                hires_bbox["ymin"]:hires_bbox["ymax"],]
aligned2 <- aligned[hires_bbox["xmin"]:hires_bbox["xmax"], 
                    hires_bbox["ymin"]:hires_bbox["ymax"],]
detected2 <- detected[hires_bbox["xmin"]:hires_bbox["xmax"], 
                      hires_bbox["ymin"]:hires_bbox["ymax"],]

writeImage(lowres2, file.path(outpath, "spatial", "tissue_lowres_image.png"))
writeImage(hires2, file.path(outpath, "spatial", "tissue_hires_image.png"))
writeImage(aligned2, file.path(outpath, "spatial", "aligned_fiducials.jpg"))
writeImage(detected2, file.path(outpath, "spatial", "detected_tissue_image.jpg"))
