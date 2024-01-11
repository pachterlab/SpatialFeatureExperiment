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


# get images ----
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


# prepare manifest.json ----
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

# get mols coords ----
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

# downsample mols a bit more
set.seed(567)
mols_joined %<>% 
  dplyr::sample_frac(size = 0.1)

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


# get count matrix ----
mat_sub <- assay(sfe_sub, "counts")
# filter count matrix to keep genes present in mols coords
gene_indx <- which(rownames(mat_sub) %in% unique(mols_joined$gene))
mat_sub <- mat_sub[gene_indx,]

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


# get metadata ----
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


# get cell segmentations ----
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


# test loading the toy dataset ----
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


# convert txSpots to hdf5 files ----
library(rhdf5)
dir_use <- "./vizgen_test_repo/vizgen_hdf5/"
# load another dataset file to use an example (local data)
fn <- file.path(dir_use, "cell_boundaries", "feature_data_119.hdf5")
hdf5_eg <- rhdf5::h5dump(fn)[[1]]

# coords of 1 cell for each z plane
hdf5_eg[[1]] %>% str
# z-plane 0
hdf5_eg[[1]]$zIndex_0$p_0$coordinates %>% str

# make a list of z planes and coords per cell
# clone example for one cell
cellsegs <- hdf5_eg[[2]]
cellsegs_out <-
  bplapply(cellSeg(sfe_sub) %>% 
             nrow() %>% seq(), function(i) {
               # get coords per cell
               coordinates <-
                 cellSeg(sfe_sub) %>%
                 dplyr::slice(i) %>%
                 st_geometry() %>%
                 st_coordinates() %>% 
                 as.data.frame() %>% 
                 # change sign to positive for y coord
                 dplyr::transmute(X = X, Y = -Y) %>% t()
               
               # convert to an array
               dim(coordinates) <- c(2, ncol(coordinates), 1)
               
               # replace coords
               cellseg <- 
                 lapply(cellsegs[-which(names(cellsegs) == "z_coordinates")], 
                        function(x) {
                          x$p_0$coordinates <- coordinates
                          return(x)}
                 )
               cellseg[["z_coordinates"]] <- cellsegs[["z_coordinates"]]
               return(cellseg)  
             }, BPPARAM = BiocParallel::MulticoreParam(10, 
                                                       tasks = 50L,
                                                       force.GC = FALSE, 
                                                       progressbar = TRUE))
names(cellsegs_out) <- cellSeg(sfe_sub) %>% rownames()
cellsegs_out[[1]] %>% str
cellsegs_out %>% length

dir_github <- "./vizgen_test_repo/vizgen_hdf5_github/"
# prepare as hdf5
dir.create(file.path(dir_github, "cell_boundaries"))
# z-planes, only 2nd and 3rd to make less larger files
z_names <- paste0("zIndex_", seq(2,3))
# in total 4 hdf5 files
new_fns <- file.path(dir_github, "cell_boundaries", 
                     paste0("feature_data_z2_z3_", seq(4), ".hdf5"))
new_fns
# use 80 cells in total, 20 per file
cell_ids <- 
  names(cellsegs_out)[seq(1, length(cellsegs_out), length.out = 80)]
cell_ids %>% str

# split into 4 hdf5 files
bplapply(seq(new_fns), function(x) {
  # make files
  h5createFile(new_fns[x])
  h5createGroup(new_fns[x], "featuredata")
  # use per cell and per zIndex
  bplapply(seq(z_names), function(z) {
    bplapply(split(seq(80), rep(seq(4), 20))[[x]], function(i) {
      gn <<- file.path("featuredata", cell_ids[i], z_names[z], "p_0/coordinates")
      h5createGroup(new_fns[x], file.path("featuredata", cell_ids[i]))
      h5createGroup(new_fns[x], file.path("featuredata", cell_ids[i], z_names[z]))
      h5createGroup(new_fns[x], file.path("featuredata", cell_ids[i], z_names[z], "p_0"))
      coords <<- cellsegs_out[[cell_ids[i]]][[z_names[z]]]$p_0$coordinates
      h5write(coords, new_fns[x], gn)
    }, BPPARAM = BiocParallel::SerialParam()) %>% invisible()
    
  }, BPPARAM = BiocParallel::SerialParam(force.GC = FALSE, 
                                         progressbar = TRUE)) %>% 
    suppressWarnings() %>% suppressMessages() %>% invisible()	
}, BPPARAM = BiocParallel::SerialParam()) %>% invisible()


# check the structure
h5ls(new_fns[1]) %>% str
test1 <- rhdf5::h5dump(new_fns[1])[[1]]
test1[[1]] %>% str
test1 %>% length

# test
polys <- SpatialFeatureExperiment:::.h52poly_fov(new_fns[1], z = 2)
polys %>% str


# test loading hdf5 files dir ----
dir_github <- "./vizgen_test_repo/vizgen_hdf5_github/"
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

# obj has 77 cells if `filter_counts = TRUE`, else 80 cells.
sfe
# normalize raw counts
sfe %<>% logNormCounts()
sfe
imgData(sfe)
rowGeometry(sfe) %>% str # XYZ coords

# plot gene expression
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
