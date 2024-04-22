## make toy objects for Seurat -> SFE converter tests

# load libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(BiocParallel)
  library(tools)
  library(rlang)
  library(data.table)
})

# set R options, paths, object params, names, functions, etc..
message("Setting options, paths and object params..")
Sys.setenv(LANG = "en")
# ..to visualize huge array, standard notation vs scientific one, digits after comma.
options(max.print = 6e+5, scipen = 500, digits = 6)

# set maximum for object size
# [Parallelization in Seurat](https://satijalab.org/seurat/articles/future_vignette.html)
# for global variables size, if default limit is smaller..
message("Setting default limit for globar vars size: ", 
        paste0("~", round((50000 * 1024^2) / 1e+9), "G RAM"))
options(future.globals.maxSize = 50000 * 1024^2)

#' To get metadata df from sfe object
#' @param object SFE or Seurat object
#' @noRd
callMeta <- function(object = NULL) {
  if (class(object) == "Seurat") {
    #is(object, "Seurat")
    # return Seurat metadata
    return(object@meta.data) 
  } else {
    SpatialFeatureExperiment::callMeta(object)
  }
}


# TODO: update input paths using data SFEData ----
# save objects in this dir
dir_out <- "./inst/extdata/"

# path to dataset
dir_use <- "./inst/extdata/vizgen_cellbound/"

# load data - Seurat -- Vizgen ----
# NOTE, offical Seurat::LoadVizgen is not compatible with the latest Vizgen outs..
#..here we use this PR to handle parquet files https://github.com/satijalab/seurat/pull/7190
#..install it using:
# Seurat v4 `remotes::install_github(repo = 'alikhuseynov/seurat', ref = 'feat/vizgen')`
# Seurat v5 `remotes::install_github(repo = 'alikhuseynov/seurat', ref = 'vizgen_seurat5')`
obj_vz <- 
  LoadVizgen(data.dir = dir_use,
             fov = "vz_toy.1",
             assay = "Vizgen",
             metadata = c("volume", "fov"), # add cell volume info
             type = c("segmentations", "centroids"), # type of cell spatial coord matrices
             add.zIndex = TRUE, # add z slice section to a cell
             update.object = TRUE,
             use.BiocParallel = TRUE,
             workers.MulticoreParam = 12, # for `BiocParallel` processing
             min.area = 5,
             add.molecules = TRUE, # add "molecules" coordinates to FOV of the object
             verbose = TRUE
  )
obj_vz
GetAssayData(obj_vz, assay = "Vizgen", slot = "data") %>% str
# for SeuratObject v5, use layer instead
GetAssayData(obj_vz, assay = "Vizgen", layer = "data") %>% str
obj_vz$orig.ident <- "vz.toy.1"

# plot it
ImageFeaturePlot(obj_vz,
                 fov = Images(obj_vz), # using cropped FOV for plotting.
                 features = rownames(obj_vz)[1],
                 molecules = rownames(obj_vz)[1],
                 mols.cols = "cyan1",
                 mols.size = 0.5,
                 #nmols = 1000,
                 size = 1,
                 border.size = NA,
                 #border.color = alpha("blue1", 0.3),
                 boundaries = "segmentation",
                 alpha = 0.7,
                 axes = TRUE,
                 combine = TRUE)

# 2nd sample would be bboxed or cropped 1st toy sample
# load modified subset function for Seurat
devtools::source_url(url = "https://github.com/alikhuseynov/add-on_R/blob/develop/R/subset_obj_seurat_v2.R?raw=TRUE")
# eg subset for 100 genes and 4 molecules
# NOTE, FOVs will be subsetted as well
set.seed(999)
random_cells <- sample(Cells(obj_vz), size = 500)
obj_vz_sub <- 
  subset_opt(obj_vz, 
             cells = random_cells)
obj_vz_sub$orig.ident <- "vz.toy.2"
obj_vz_sub
# change FOV name
names(obj_vz_sub@images) <- "vz.toy.2"

# subset the 1st object to remove cells of 2nd obj
obj_vz_1 <- 
  subset_opt(obj_vz, 
             cells = which(!Cells(obj_vz) %in% random_cells))
obj_vz_1$orig.ident <- "vz.toy.1"
obj_vz_1

# merge 2 objects
obj_vz_merged <- merge(obj_vz_1, obj_vz_sub)
obj_vz_merged
# export obj - 1 sample
saveRDS(obj_vz, file.path(dir_out, "seu_vz_toy.rds"))
# ..2 samples
saveRDS(obj_vz_merged, file.path(dir_out, "seu_vz_toy_multi.rds"))


# load data - Seurat -- Xenium ----
# same data as shared full `xenium_toy`..
#..Xenium standard test (in-house) FFPE, Jurkat and Raji cell mix - run Nov 2023
# custom path to dataset
dir_use <- "./test/test_xenium_run/"

# Seurat doesn't load cell_feature_matrix.h5
#.. this might be fixed in recent PRs (not merged yet)
# alternative, load data to get cell ids etc..
xen_toy_data <- 
  ReadXenium(data.dir = dir_use, 
             outs = "microns",
             type = c("centroids", "segmentations"))
xen_toy_data %>% str

# load object
obj_xen <- 
  LoadXenium(data.dir = dir_use, 
             assay = "Xenium",
             fov = "xen.toy.1")
obj_xen
# subset to keep cells from toy dataset
obj_xen <- 
  subset_opt(obj_xen, 
             cells = xen_toy_data[[6]])
obj_xen
GetAssayData(obj_xen, assay = "Xenium", slot = "data") %>% str
# for SeuratObject v5, use layer instead
GetAssayData(obj_xen, assay = "Xenium", layer = "data") %>% str
obj_xen$orig.ident <- "xen.toy.1"

# plot it
ImageFeaturePlot(obj_xen,
                 fov = Images(obj_xen), # using cropped FOV for plotting.
                 features = rownames(obj_xen)[1],
                 molecules = rownames(obj_xen)[1],
                 mols.cols = "cyan1",
                 mols.size = 0.5,
                 #nmols = 1000,
                 size = 1,
                 border.size = NA,
                 #border.color = alpha("blue1", 0.3),
                 boundaries = "segmentation",
                 alpha = 0.7,
                 axes = TRUE,
                 combine = TRUE)

# 2nd sample would be bboxed or cropped 1st toy sample
# NOTE, FOVs will be subsetted as well
set.seed(999)
random_cells <- sample(Cells(obj_xen), size = 5000)
obj_xen_sub <- 
  subset_opt(obj_xen, 
             cells = random_cells)
obj_xen_sub$orig.ident <- "xen.toy.2"
obj_xen_sub
# change FOV name
names(obj_xen_sub@images) <- "xen.toy.2"

# subset the 1st object to remove cells of 2nd obj
obj_xen_1 <- 
  subset_opt(obj_xen, 
             cells = which(!Cells(obj_xen) %in% random_cells))
obj_xen_1$orig.ident <- "xen.toy.1"
obj_xen_1

# merge 2 objects
obj_xen_merged <- merge(obj_xen_1, obj_xen_sub)
obj_xen_merged
# export obj - 1 sample
saveRDS(obj_xen, file.path(dir_out, "seu_xen_toy.rds"))
# ..2 samples
saveRDS(obj_xen_merged, file.path(dir_out, "seu_xen_toy_multi.rds"))


# load data - Seurat -- Visium ----
# path to dataset
dir_use <- "./inst/extdata/"
#samples: `sample01` & `sample02`
files_vis <- 
  list.dirs(dir_use, full.names = TRUE) %>% 
  grep("outs$", ., value = TRUE)

# Since Seurat doesn't support reading MTX files for Visium
# https://github.com/satijalab/seurat/issues/5806
# Load in expression data and create a Seurat object
obj_vis_list <-
  lapply(seq(files_vis), function(i) {
    counts <- Read10X(file.path(files_vis[i], "filtered_feature_bc_matrix"))
    obj_vis <- CreateSeuratObject(counts)
    # Load in the image data
    img <- Read10X_Image(file.path(files_vis[i], "spatial"))
    # Correct the image data to match the Seurat object
    img <- img[Cells(obj_vis)]
    DefaultAssay(img) <- DefaultAssay(obj_vis)
    # Add the image to the object
    obj_vis[[paste0("sample0", i)]] <- img
    return(obj_vis)
  })
obj_vis_list
# quick plot
SpatialFeaturePlot(obj_vis_list[[1]], features = "nCount_RNA")

# merge 2 objects
obj_vis <- merge(obj_vis_list[[1]], obj_vis_list[[2]])
obj_vis
# export obj - 1 sample
saveRDS(obj_vis_list[[1]], file.path(dir_out, "seu_vis_toy.rds"))
# ..2 samples
saveRDS(obj_vis, file.path(dir_out, "seu_vis_toy_multi.rds"))


# TODO load data - Seurat -- Visium HD ----
# NOTE: 
# one needs SeuratObject "5.0.1.9008"
#remotes::install_github(repo = "satijalab/seurat-object", ref = "develop")
# Seurat "visium-hd" branch "5.0.3.9909"
#remotes::install_github(repo = "satijalab/seurat", ref = "visium-hd")
# path to Visium HD directory
dir_use <- "~/Downloads/Visium_HD_Mouse_Brain/" 
# it must contain these dirs:
#../Visium_HD_Mouse_Brain
#  └── binned_outputs
#  ├── square_008um
#  └── square_016um
# see this for details 
# https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview#hd-outputs
# dataset used is -> https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he

# example is -> https://satijalab.org/seurat/articles/visiumhd_commands_intro#load-visium-hd-data
# loading 2 spatially binned resolutions 8 and 16um
obj_hd <- Load10X_Spatial(data.dir = dir_use, bin.size = c(8, 16))
# results in 2 assays and 2 samples
obj_hd
# plot counts 
SpatialFeaturePlot(obj_hd, features = "nCount_Spatial.008um", pt.size.factor = 4)

# subset object to keep only first 50 genes
obj_hd %<>% subset(features = rownames(obj_hd)[1:50])
DefaultAssay(obj_hd) <- "Spatial.008um"
# log Normalize data to set "data" layer or slot
obj_hd %<>% NormalizeData()

# save object with 2 FOVs
saveRDS(obj_hd, file.path(dir_out, "seu_vishd_multi.rds"))

