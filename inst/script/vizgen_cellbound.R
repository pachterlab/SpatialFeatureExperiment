#
# load libs - less libs
suppressPackageStartupMessages({
  library(ggplot2)
  #library(Seurat)
  library(dplyr)
  library(magrittr)
  library(BiocParallel)
  library(scuttle)
  #library(SingleCellExperiment) 
  #library(SpatialExperiment) 
  library(SpatialFeatureExperiment)
  library(Voyager)
  library(terra)			
  library(sf)
})

## ------------------------- ##
## Brain cancer toy dataset 1
## ------------------------- ##
# my dir to large test dataset
dir_use <- "./vizgen_test_repo/vizgen_cellbound/"
dir_github <- "./vizgen_test_repo/vizgen_cellbound_github/"

# load large SFE object
sfe <-
  readVizgen(data_dir = dir_use,
             z = "all",
             z_option = "3d", # this will return XYZ MULTIPOINT for rowGeometiries
             sample_id = "vizgen_test",
             min_area = 15,
             image = c("DAPI", "PolyT", "Cellbound"),
             flip = "geometry", # "image" & "none", "geometry"
             max_flip = "50 MB",
             filter_counts = TRUE, # keep cells with counts > 0 or not    
             add_molecules = TRUE,
             use_bboxes = FALSE,
             #file_out = file.path(dir_use, "detected_transcripts.parquet"),
             BPPARAM = BiocParallel::MulticoreParam(14, 
                                                    tasks = 80L,
                                                    force.GC = FALSE, 
                                                    progressbar = TRUE)
  )
sfe
# normalize raw counts
sfe %<>% logNormCounts()
sfe
imgData(sfe)
rowGeometry(sfe) %>% str # XYZ coords

# check if flip works
colGeometry(sfe, 1) %>% st_geometry() %>% st_bbox
cellSeg(sfe) %>% st_geometry() %>% st_bbox
txSpots(sfe) %>% st_geometry() %>% st_bbox

# plot it
# Segs
options(repr.plot.height = 5, repr.plot.width = 10)
pl1 <- 
  plotSpatialFeature(sfe, features = "COL1A2", 
                     size = 4,
                     #colGeometryName = "centroids", 
                     colGeometryName = "cellSeg",
                     dark = TRUE,
                     #scattermore = TRUE, # will plot only centroids!
                     image_id = "Cellbound2_z3" # "DAPI_z3"
  ) & Seurat::DarkTheme()
pl1

# Segs using bbox cropping
bbox_use <- c(xmin = 6500, ymin = -1500, xmax = 6800, ymax = -1200)
plotSpatialFeature(sfe, features = "COL1A2", bbox = bbox_use,
                   size = 4, 
                   #colGeometryName = "centroids",
                   colGeometryName = "cellSeg",
                   dark = TRUE,
                   #scattermore = TRUE, # will plot only centroids!
                   image_id = c("Cellbound2_z3"),
) & Seurat::DarkTheme()
#theme(axis.line = element_line(color = "white"), 
#      axis.ticks = element_line(color = "white", size = 1, linetype = 1),
#      axis.text = element_text(color = "white"))

# subset obj given the bbox
sfe_sub <- 
  SpatialFeatureExperiment::crop(sfe, 
                                 colGeometryName = "cellSeg",
                                 sample_id = "vizgen_test",
                                 y = bbox_use)
# keep only background genes and some panel general markers
sfe_sub <- sfe_sub[c(1:50, grep("Blank-", rownames(sfe_sub)))]
sfe_sub

# plot it - Segs
options(repr.plot.height = 5, repr.plot.width = 10)
pl1 <- 
  plotSpatialFeature(sfe_sub, features = "COL1A2", 
                     size = 4,
                     #colGeometryName = "centroids", 
                     colGeometryName = "cellSeg",
                     dark = TRUE,
                     #scattermore = TRUE, # will plot only centroids!
                     image_id = "Cellbound2_z3" # "DAPI_z3"
  ) & Seurat::DarkTheme()
pl1

# check if flip stays
colGeometry(sfe_sub, 1) %>% st_geometry() %>% st_bbox
cellSeg(sfe_sub) %>% st_geometry() %>% st_bbox
txSpots(sfe_sub) %>% st_geometry() %>% st_bbox

# export sfe_toy obj
dir_local <- "./SpatialFeatureExperiment/seurat_v4/inst/extdata"
saveRDS(sfe_sub, file = file.path(dir_local, "sfe_vizgen_toy.rds"))
#sfe_sub <- readRDS(file = file.path(dir_local, "sfe_vizgen_toy.rds"))
sfe_sub

# OK -> get images ----
imgData(sfe_sub)
imgData(sfe_sub)[1,] # 1st image
# get images in a list
im_sub <- 
  lapply(imgData(sfe_sub) %>% nrow() %>% seq(), 
         function(i) {
           im_sub <- 
             Voyager:::.get_img_df(sfe_sub, 
                                   sample_id = "vizgen_test", 
                                   image_id = imgData(sfe_sub)$image_id[i], 
                                   bbox = NULL # set bbox arg to NULL
             )
         })
names(im_sub) <- imgData(sfe_sub)$image_id

# get image - in a loop
img <- 
  lapply(seq(im_sub), function (i) {
    im_sub[[i]]$data[[1]]@image %>% 
      unwrap() # convert from PackedSpatRaster to SpatRaster
  })
plot(img[[1]])
# convert to export as image
img_out.list <-
  lapply(seq(img), function(i) {
    # using terra lib
    img_out <- 
      img[[i]] %>%
      terra::aggregate(., fact = 4) %>%
      #terra::flip() %>%
      terra::as.array()
    img_out[is.nan(img_out)] <- 0
    return(img_out)
  })
# export images
library(tiff)
for (i in seq(img_out.list)) {
  writeTIFF(img_out.list[[i]] / c(max(img_out.list[[i]]) * 3), # downsample image
            paste0(dir_github, "images/mosaic_", names(im_sub)[i], ".tif"))
}
# load image and plot it
img_test <- rast(file.path(dir_github, "images/mosaic_Cellbound1_z3.tif"))
plot(img_test)

# OK -> prepare manifest.json ----
library(jsonlite)
manifest <- read_json(file.path(dir_use, "images/manifest.json"), simplifyVector = TRUE)
manifest %>% str
mnfst <- manifest
# keep only image names present in downsampled obj
mnfst$mosaic_files %<>% 
  filter(file_name %in% paste0("mosaic_", names(im_sub), ".tif"))
mnfst$mosaic_pyramid_files <- NULL
# extract relevant image dims
img_out <- 
  img[[1]] %>% 
  terra::aggregate(., fact = 4) #%>% terra::flip()
mnfst$mosaic_width_pixels <- ncol(img_out)
mnfst$mosaic_height_pixels <- nrow(img_out)
# make bbox with positive signes
extent <- as.vector(ext(img_out))[c("xmin", "ymin", "xmax", "ymax")] 
extent[c("ymin", "ymax")] <- -extent[c("ymax", "ymin")]	
extent %>% unname()
mnfst$bbox_microns <- extent %>% unname()
mnfst$hor_num_tiles_box <- 1
mnfst$vert_num_tiles_box <- 1
mnfst %>% str
write_json(mnfst, file.path(dir_github, "images/manifest.json"), 
           pretty = TRUE,
           auto_unbox = TRUE)

# OK -> get count matrix ----
mat_sub <- assay(sfe_sub, "counts")
# read original count matrix
mat <- data.table::fread(file.path(dir_use, "cell_by_gene.csv"), 
                         colClasses = list(character = 1))
mat %>% str
# match them
mat %<>%
  # keep cells from mat_sub
  dplyr::slice(match(mat_sub %>% colnames(), mat$cell)) %>%
  # keep genes from mat_sub, and cell col
  select(c(cell, 
           match(mat_sub %>% rownames(), mat %>% colnames())))
# export it
data.table::fwrite(mat, file = file.path(dir_github, "cell_by_gene.csv"))

# OK -> get metadata ----
meta_orig <- data.table::fread(file.path(dir_use, "cell_metadata.csv"))
meta_orig %>% str

# get metadata df from sce, sfe or miloR obj
callMeta <- function(object = NULL) {
  return(colData(object)@listData %>% as.data.frame.list())
}
meta_sub <- 
  callMeta(sfe_sub) %>% 
  select(-contains("sample"))
meta_sub %>% str

# subset original metadata given cell ids
cells_use <- 
  match(meta_sub %>% rownames(), meta_orig$EntityID %>% as.character())
# check if cell ids correspond
identical(meta_orig %>%
            dplyr::slice(cells_use) %>% 
            pull(EntityID) %>% as.character, 
          meta_sub %>% rownames())
# all good!
meta_sub <-
  meta_orig %>%
  dplyr::slice(cells_use)
# export it
data.table::fwrite(meta_sub, file = file.path(dir_github, "cell_metadata.csv"))

# OK -> get cell segmentations ----
parq_orig <- sfarrow::st_read_parquet(file.path(dir_use, "cell_boundaries.parquet"))
parq_orig %>% str
cellSeg(sfe_sub) %>% str

# flip y coords
mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
parq <- cellSeg(sfe_sub)
st_geometry(parq) <- c(st_geometry(parq) * mat_flip)
st_geometry(parq) %>% st_bbox

# Not using that field anyway so can be random
parq$ID <- sample(seq_len(nrow(parq)) - 1, nrow(parq))
parq$EntityID <- 
  cellSeg(sfe_sub) %>% 
  rownames() %>% 
  bit64::as.integer64.character(.)
#parq$ZIndex <- 0
#parq$Type <- "cell"
#parq$ZLevel <- 1.5
parq$ParentID <- parq$ParentType <- parq$Name <- NA
parq$X__index_level_0__ <- parq$ID
names(parq)[names(parq) == "geometry"] <- "Geometry"
st_geometry(parq) <- "Geometry"
parq <- parq[,names(parq_orig), drop = FALSE]
parq %>% str
# export file
sfarrow::st_write_parquet(parq, 
                          file.path(dir_github, "cell_boundaries.parquet")) %>% suppressWarnings()

# OK -> get mols coords ----
# XYZ coords
rowGeometry(sfe_sub) %>% str
rowGeometry(sfe_sub) %>% st_geometry() %>% str
# flip mols coords
txSpots(sfe_sub) %>% st_geometry() %>% st_bbox()
mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
(st_geometry(txSpots(sfe_sub)) * mat_flip) %>% st_bbox
mols <- txSpots(sfe_sub)
st_geometry(mols) <- (st_geometry(mols) * mat_flip)
# crop before convering mols
mols <-
  st_crop(mols, 
          cellSeg(sfe_sub) %>%
            st_geometry() %>%
            st_bbox())
mols %>% str

# convert mols
mols <-
  bplapply(rownames(sfe_sub) %>% seq(), function(i) {
    mols %>%
      dplyr::slice(i) %>%
      st_geometry() %>%
      st_coordinates() %>%
      as.data.frame() %>%
      dplyr::transmute(global_x = X, 
                       global_y = Y,
                       #global_z = Z,
                       gene = rownames(sfe_sub)[i])
  }, BPPARAM = BiocParallel::MulticoreParam(12, 
                                            tasks = 50L,
                                            force.GC = FALSE, 
                                            progressbar = TRUE)
  ) %>% do.call(bind_rows, .)
mols %>% str

# load original molecule coords
mols_orig <- data.table::fread(file.path(dir_use, "detected_transcripts.csv"))
mols_orig %>% str

# plot 1 molecule
mols %>%
  filter(gene == "COL1A2") %>%
  ggplot(aes(global_x, global_y)) &
  geom_hex(bins = 50)

# filter orginal mols given cropped mols range
bbox_mols <-
  mols %>% 
  select(contains("global")) %>%
  apply(., 2, range)
bbox_mols
# filter
mols_filt <-
  filter(mols_orig %>% mutate(global_y = -global_y),
         between(global_x, bbox_mols[1], bbox_mols[2]) & 
           between(global_y, bbox_mols[3], bbox_mols[4]))
mols_filt %>% str
# join dfs
mols_joined <-
  dplyr::inner_join(mols, mols_filt)
mols_joined %>% str
# check range
mols_joined %>% 
  select(contains("global")) %>%
  apply(., 2, range)

# plot 1 molecule for all z-planes
mols_joined %>%
  filter(gene == "COL1A2") %>%
  ggplot(aes(global_x, global_y, color = global_z)) &
  geom_point(shape = 3) & 
  # add previous plot
  pl1
# things seem to correspond, woo!

# order vars
mols_joined %<>%
  # make y coord positive
  mutate(global_y = -global_y) %>%
  select(., names(mols_filt))
mols_joined %>% str

# export transcripts coords
data.table::fwrite(mols_joined, file.path(dir_github, "detected_transcripts.csv"))

# OK -> test loading the toy dataset ----
dir_github <- "./vizgen_test_repo/vizgen_cellbound_github/"
# load SFE object
sfe <-
  readVizgen(data_dir = dir_github,
             z = "all",
             z_option = "3d", # this will return XYZ MULTIPOINT for rowGeometiries
             sample_id = "vizgen_toy",
             min_area = 15,
             image = c("DAPI", "PolyT", "Cellbound"),
             flip = "geometry", # "image" & "none", "geometry"
             max_flip = "50 MB",
             filter_counts = TRUE, # keep cells with counts > 0 or not    
             add_molecules = TRUE,
             use_bboxes = FALSE,
             #file_out = file.path(dir_use, "detected_transcripts.parquet"),
             BPPARAM = BiocParallel::MulticoreParam(10, 
                                                    tasks = 50L,
                                                    force.GC = FALSE, 
                                                    progressbar = TRUE)
  )
sfe
# normalize raw counts
sfe %<>% logNormCounts()
sfe
imgData(sfe)
rowGeometry(sfe) %>% str # XYZ coords

# toy obj
colGeometry(sfe, 1) %>% st_geometry() %>% st_bbox
cellSeg(sfe) %>% st_geometry() %>% st_bbox
txSpots(sfe) %>% st_geometry() %>% st_bbox

# plot it
# Segs
options(repr.plot.height = 5, repr.plot.width = 10)
pl1 <- 
  plotSpatialFeature(sfe, features = "COL1A2", 
                     size = 4,
                     #colGeometryName = "centroids", 
                     colGeometryName = "cellSeg",
                     dark = TRUE,
                     #scattermore = TRUE, # will plot only centroids!
                     image_id = "Cellbound2_z3", #"DAPI_z3",
  ) & Seurat::DarkTheme()
pl1
