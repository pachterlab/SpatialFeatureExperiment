## make Xenium toy dataset for `readXenium` ----

# load libs - less libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  #library(Seurat)
  library(dplyr)
  library(magrittr)
  library(BiocParallel)
  library(scuttle)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  #library(SpatialFeatureExperiment)
  library(Voyager)
  library(terra)
  library(sf)
  library(tools)
  library(rlang)
  library(data.table)
  library(DropletUtils)
  library(RBioFormats)
  library(jsonlite)
})
devtools::load_all()
## ------------------------------------- ##
## Xenium standard test (in-house) FFPE
## Jurkat and Raji cell mix - run Nov 2023
## ------------------------------------- ##
# laod xenium test data ----
data_dir <- "xenium_toy"
# out dir
dir_github <- "xenium_github"
#BPPARAM <- BiocParallel::MulticoreParam(8, tasks = 50L,
#                                        force.GC = FALSE,
#                                        progressbar = TRUE)
BPPARAM <- MulticoreParam(3, progressbar = TRUE)

sfe <-
  readXenium(data_dir = data_dir,
             sample_id = "test_xenium",
             image = c("morphology_focus", "morphology_mip"),
             segmentations = c("cell", "nucleus"),
             row.names = "symbol",
             read.image_args = # list of arguments to passed to RBioFormats::read.image
               list("resolution" = 4L,
                    "filter.metadata" = TRUE,
                    "read.metadata" = FALSE,
                    "normalize" = FALSE),
             image_threshold = 30,
             flip = "geometry",
             filter_counts = FALSE,
             add_molecules = TRUE,
             BPPARAM = BPPARAM,
             file_out = file.path(data_dir, "tx_spots.parquet")
  )
# Cells seem pretty uniformly distributed, not surprising since it's cell culture
plotGeometry(sfe, "cellSeg")
plotCellBin2D(sfe, bins = 100, hex = TRUE)
sfe$nCounts <- colSums(counts(sfe))
plotSpatialFeature(sfe, "nCounts", colGeometryName = "cellSeg")
plotSpatialFeature(sfe, "nCounts", colGeometryName = "nucSeg")

# normalize raw counts
#sfe %<>% logNormCounts()
#sfe
#imgData(sfe)
#rowGeometry(sfe) %>% str

# Segs using bbox cropping
bbox_use <- c(xmin = 1000, ymin = -1400, xmax = 1500, ymax = -1000)
plotSpatialFeature(sfe,
                   features = "nCounts",
                   #features = "NegControlCodeword_0500",
                   bbox = bbox_use,
                   #size = 4,
                   #colGeometryName = "centroids",
                   colGeometryName = "nucSeg"#,
                   #dark = TRUE,
                   #scattermore = TRUE, # will plot only centroids!
                   #image_id = imgData(sfe)$image_id[1],
) #& Seurat::DarkTheme()

# make a small cropped object ----
bbox_use <- c(xmin = 1000, ymin = -1400, xmax = 1400, ymax = -1000)
sfe_sub <-
  SpatialFeatureExperiment::crop(sfe,
                                 colGeometryName = "nucSeg",
                                 sample_id = "test_xenium",
                                 y = bbox_use)
sfe_sub
pl1 <-
  plotSpatialFeature(sfe_sub,
                     features = rownames(sfe_sub)[1],
                     #features = "NegControlCodeword_0500",
                     #bbox = bbox_use,
                     size = 4,
                     #colGeometryName = "centroids",
                     colGeometryName = "nucSeg",
                     dark = TRUE,
                     #scattermore = TRUE, # will plot only centroids!
                     image_id = imgData(sfe_sub)$image_id[1],
  ) #& Seurat::DarkTheme()
pl1


# get images ----
# read image metadata
img_fns <- list.files(data_dir,
                      full.names = TRUE,
                      pattern = "ome.tif") %>%
  grep("morphology_", ., value = TRUE)
img_fns
meta_im <-
  lapply(seq(img_fns), function(i) {
    read.omexml(img_fns[i]) %>%
      XML::xmlInternalTreeParse() %>%
      XML::xmlToList()
    #strsplit(., split = " ") %>% unlist %>%
    #grep("PhysicalSize|Size[XY]", ., value = TRUE)
  })
# image metadata
meta_im %>% str

# read the actual image with low resolution
img_fns <- list.files(data_dir,
                      full.names = TRUE,
                      pattern = "ome.tif") %>%
  grep("morphology_", ., value = TRUE)
img_fns
imgs <-
  lapply(seq(img_fns), function(i)
    read.image(img_fns[i],
               #subset = list(x = 1200:1300, y = 1200:1300),
               resolution = 4L,
               #proprietary.metadata = TRUE,
               filter.metadata = TRUE,
               read.metadata = TRUE,
               normalize = FALSE))
imgs %>% str
# plot one image
imgs[[1]]@.Data %>% t %>% rast %>% plot

# from object, get images in a list
imgData(sfe_sub)
im_sub <-
  lapply(imgData(sfe_sub) %>% nrow() %>% seq(),
         function(i) {
           im_sub <-
             Voyager:::.get_img_df(sfe_sub,
                                   sample_id = "test_xenium",
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
img[[1]] %>% plot

# Add bbox-adjested full image dims!
# full image extent
#c(0, full_res_img_width * px_size_micron,
#  0, full_res_img_heigth * px_size_micron)
# $ PhysicalSizeY    : chr "0.2125"
# $ PhysicalSizeX    : chr "0.2125"
# $ SizeY            : chr "10307"
# $ SizeX            : chr "11454"
extent <-
  setNames(c(0, 11454 * 0.2125,
             0, 10307 * 0.2125),
           c("xmin", "xmax", "ymin", "ymax"))
img_test <- imgs[[1]]@.Data %>% t %>% rast
ext(img_test) <- extent
# crop using bbox
img_crop <-
  terra::crop(terra::flip(img_test),
              c("xmin" = 1000, "xmax" = 1400, "ymin" = 1000, "ymax" = 1400))
img_crop <- terra::flip(img_crop)
img_crop %>% plot
img_crop
extent <- ext(img_crop) %>% as.vector
# flip it
extent[c("ymin", "ymax")] <- -extent[c("ymax", "ymin")]
ext(img_crop) <- extent
c(extent / 0.2125) # cropped pixel image dim
#xmin: 4706.46540880503 xmax: 6587.45073375262 ymin: -6585.91692546584 ymax: -4705.36956521739

# adjust image given the cropped object bbox
imgs %>% str
for (i in seq(imgs)) {
  imgs[[i]]@.Data <-
    matrix(values(img[[i]]) %>% as.numeric,
           ncol = ncol(img[[i]]),
           nrow = nrow(img[[i]]))
  coreMetadata(imgs[[i]])$sizeX <- ncol(img[[i]])
  coreMetadata(imgs[[i]])$sizeY <- nrow(img[[i]])
  # add part of original metadata
  globalMetadata(imgs[[i]]) <- meta_im[[i]]$Image$Pixels$.attrs
  # change full image dim to low res cropped image
  globalMetadata(imgs[[i]])[c("SizeX", "SizeY")] <-
    sqrt(c(extent / 0.2125)[c("xmax", "ymax")] ^ 2)
}

# exporting as `.ome.tif` with single res. ----
for (i in seq(imgs)) {
  write.image(x = imgs[[i]], force = TRUE,
              file = paste0(dir_github, basename(img_fns[i])))
}
# read file to test
img_test <- read.image(paste0(dir_github, basename(img_fns[1])),
                       resolution = 1, # NOTE, only single res. is available
                       filter.metadata = TRUE,
                       read.metadata = TRUE,
                       normalize = FALSE)
img_test %>% plot


# get mols coords ----
# read original transcripts
fn_mols <- SpatialFeatureExperiment:::.check_xenium_fns(data_dir, "transcripts")
mols_parq <- arrow::read_parquet(gsub(".csv.gz", ".parquet", fn_mols))
mols_parq %>% head
mols_orig <- fread(fn_mols)
mols_orig %>% head

# transcrip coords from cropped obj
# NOTE - molecules are not cropped (as in full obj)
rowGeometry(sfe_sub) %>% str
rowGeometry(sfe_sub) %>% st_geometry() %>% str
# chack xy range
mols_orig %>%
  select(starts_with(c("x_loc", "y_loc"))) %>%
  apply(., 2, range)
# flip mols coords to match cell segs
rowGeometry(sfe_sub, 1) %>% st_geometry() %>% st_bbox()
mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
(st_geometry(rowGeometry(sfe_sub, 1)) * mat_flip) %>% st_bbox
mols <- rowGeometry(sfe_sub, 1)
st_geometry(mols) <- (st_geometry(mols) * mat_flip)
mols %>% str
# remove any NAs
mols %<>% na.omit

# plot 1 molecule
mols_orig %>%
  filter(feature_name == "CD3D") %>%
  ggplot(aes(x_location, y_location)) &
  geom_hex(bins = 50)

# filter orginal mols given cell segs bbox
bbox_mols <-
  cellSeg(sfe_sub) %>%
  st_geometry() %>%
  st_bbox()
# filter molecules
mols_filt <-
  filter(mols_orig %>% mutate(y_location = -y_location),
         between(x_location, bbox_mols[1], bbox_mols[3]) &
           between(y_location, bbox_mols[2], bbox_mols[4]) &
           # qv/min_phred >= 20
           qv >= 20)
mols_filt %>% str
mols_filt$qv %>% range

# downsample mols a bit more
set.seed(567)
mols_filt %<>%
  dplyr::sample_frac(size = 0.12) %>%
  filter(# keep genes that are present in rowGeometry
    feature_name %in% mols$ID)
# check range
mols_filt %>%
  select(contains("location")) %>%
  apply(., 2, range)
mols_filt %>% str

# plot 1 molecule for all z-planes
mols_filt %>%
  filter(feature_name == "CD3D") %>%
  ggplot(aes(x_location, y_location, color = z_location)) &
  geom_point(shape = 19) &
  # add previous plot
  pl1
# things seem to correspond, woooh!

# export transcripts coords
mols_filt %<>%
  # flip y coord back
  mutate(y_location = -y_location)
fwrite(mols_filt, file.path(dir_github, "transcripts.csv.gz"))
arrow::write_parquet(mols_filt, file.path(dir_github, "transcripts.parquet"))


# get count matrix ----
mat_sub <- assay(sfe_sub, "counts")
# filter count matrix to keep genes present in mols coords
gene_indx <- which(rownames(mat_sub) %in% unique(mols_filt$feature_name))
mat_sub <- mat_sub[gene_indx,]

# export it
DropletUtils::write10xCounts(path = file.path(dir_github, "cell_feature_matrix.h5"),
                             mat_sub,
                             gene.id = rowData(sfe_sub)$ID[gene_indx],
                             gene.symbol = rowData(sfe_sub)$Symbol[gene_indx],
                             gene.type = rowData(sfe_sub)$Type[gene_indx],
                             overwrite = TRUE,
                             version = "3"
)

# check it
sce <- read10xCounts(file.path(dir_github, "cell_feature_matrix.h5"),
  col.names = TRUE, row.names = c("symbol"))
identical(mat_sub %>% rownames,
          assay(sce, "counts") %>% rownames) # TRUE


# get metadata ----
meta_orig <- fread(file.path(data_dir, "cells.csv.gz"))
meta_orig %>% str
# get metadata df from sce, sfe or miloR obj
meta_sub <-
  callMeta(sfe_sub) %>%
  select(-contains("sample"))
meta_sub %>% str

# subset original metadata given cell ids
cells_use <-
  match(meta_sub %>% rownames(), meta_orig$cell_id %>% as.character())
# check if cell ids correspond
identical(meta_orig %>%
            dplyr::slice(cells_use) %>%
            pull(cell_id) %>% as.character,
          meta_sub %>% rownames())
# all good, woooh!
meta_sub <-
  meta_orig %>%
  dplyr::slice(cells_use)
# export it
data.table::fwrite(meta_sub, file = file.path(dir_github, "cells.csv.gz"))
arrow::write_parquet(meta_sub, file.path(dir_github, "cells.parquet"))
system(paste0("ls -lth ", dir_github), intern = TRUE)


# get cell segmentations ----
fn_segs <- list.files(data_dir, pattern = "boundaries.", full.names = TRUE)
fn_segs
polys_orig <- lapply(grep("boundaries.csv.gz", fn_segs, value = TRUE), fread)
polys_orig %>% str
parq_orig <- lapply(grep("boundaries.parquet", fn_segs, value = TRUE),
                    arrow::read_parquet)
parq_orig %>% head

# subset segs using object's cells
polys_sub <-
  lapply(seq(polys_orig), function(i) {
    polys_orig[[i]] %>%
      filter(cell_id %in% colnames(sfe_sub))
  })
polys_sub %>% str
identical(polys_sub[[1]]$cell_id %>% unique, colnames(sfe_sub))

# subset parquet given cell indices in polys_orig
parq_sub <-
  lapply(seq(parq_orig), function(i) {
    cell_index <-
      which(polys_orig[[i]]$cell_id %in% colnames(sfe_sub))
    parq_orig[[i]] %>%
      dplyr::slice(cell_index)
  })
parq_sub %>% str

# export .csv.gz & .parquet
for (i in seq(polys_sub)) {
  fwrite(polys_sub[[i]], file = file.path(dir_github,
                                          grep("boundaries.csv.gz", basename(fn_segs), value = TRUE))[i])
  arrow::write_parquet(parq_sub[[i]],
                       file.path(dir_github,
                                 grep("boundaries.parquet", basename(fn_segs), value = TRUE))[i])
}

# Test loading the toy dataset ----
dir_github <- "./test/xenium_github/"
# load SFE object
sfe <-
  readXenium(data_dir = dir_github,
             sample_id = "test_xenium",
             image = c("morphology_focus", "morphology_mip"),
             segmentations = c("cell", "nucleus"),
             row.names = "symbol",
             read.image_args = # list of arguments to passed to RBioFormats::read.image
               list("resolution" = 1L,
                    "filter.metadata" = TRUE,
                    "read.metadata" = FALSE,
                    "normalize" = FALSE),
             image_threshold = 30,
             flip = "geometry",
             filter_counts = TRUE,
             add_molecules = TRUE,
             BPPARAM = BPPARAM,
             file_out = NULL
             #file_out = file.path(data_dir, "tx_spots.parquet")
  )
# normalize raw counts
sfe %<>% logNormCounts()
sfe
imgData(sfe)
rowGeometry(sfe) %>% str

# compare objs
colGeometry(sfe_sub, 1) %>% st_geometry() %>% st_bbox
cellSeg(sfe_sub) %>% st_geometry() %>% st_bbox
txSpots(sfe_sub) %>% st_geometry() %>% st_bbox

# toy obj
colGeometry(sfe, 1) %>% st_geometry() %>% st_bbox
cellSeg(sfe) %>% st_geometry() %>% st_bbox
txSpots(sfe) %>% st_geometry() %>% st_bbox

# these plot must correspond! ----
# plot it
# Segs
options(repr.plot.height = 5, repr.plot.width = 10)
plotSpatialFeature(sfe,
                   features = rownames(sfe)[1],
                   #features = "NegControlCodeword_0500",
                   #bbox = bbox_use,
                   size = 4,
                   #colGeometryName = "centroids",
                   colGeometryName = "nucSeg",
                   dark = TRUE,
                   #scattermore = TRUE, # will plot only centroids!
                   image_id = imgData(sfe)$image_id[1],
) & Seurat::DarkTheme()

# downsampled obj
plotSpatialFeature(sfe_sub,
                   features = rownames(sfe_sub)[1],
                   #features = "NegControlCodeword_0500",
                   #bbox = bbox_use,
                   size = 4,
                   #colGeometryName = "centroids",
                   colGeometryName = "nucSeg",
                   dark = TRUE,
                   #scattermore = TRUE, # will plot only centroids!
                   image_id = imgData(sfe_sub)$image_id[1],
) & Seurat::DarkTheme()
