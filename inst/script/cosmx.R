library(data.table)
library(Matrix)
library(tidyverse)
library(sf)

# Experiment----------
spots <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_tx_file.csv")
unique(spots$z)
meta <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_metadata_file.csv")
polys <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter-polygons.csv")
# I think 0 means not in any cell
mat <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_exprMat_file.csv")

# Looks like cell ID's are only unique within FOV.
polys |> 
    filter(cellID == 1L, fov = 1) |> 
    ggplot(aes(x = x_global_px, y = y_global_px)) +
    geom_path() + coord_equal()

length(unique(polys$cellID)) # Less than total number of cells expected from metadata

setdiff(mat$cell_ID, meta$cell_ID)
setdiff(mat$cell_ID, unique(polys$cellID))
mat0 <- mat[mat$cell_ID == 0,] # Has row for spots outside cells in each FOV separately
# What to do with those? What are they used for? I don't think anyone would include
# them in analyses involving cells. So maybe rowData?
meta$cell_ID <- paste(meta$cell_ID, meta$fov, sep = "_")
mat$cell_ID <- paste(mat$cell_ID, mat$fov, sep = "_")
polys$cellID <- paste(polys$cellID, polys$fov, sep = "_")

mat <- mat[match(meta$cell_ID, mat$cell_ID),]
cell_ids <- mat$cell_ID
mat <- mat[,3:ncol(mat)] |> 
    as.matrix() |> 
    as("CsparseMatrix") |> t()
colnames(mat) <- cell_ids
polys <- df2sf(polys, spatialCoordsNames = c("x_global_px", "y_global_px"),
               geometryType = "POLYGON",
               id_col = "cellID")
areas <- st_area(polys)
summary(areas)
polys <- polys[match(meta$cell_ID, polys$ID),]
plot(polys$geometry[100])
all(st_is_valid(polys))

sfe <- SpatialFeatureExperiment(list(counts = mat), colData = meta,
                                spatialCoordsNames = c("CenterX_global_px", "CenterY_global_px"))
cellSeg(sfe) <- polys

spots_sf <- formatTxSpots("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_tx_file.csv",
                       spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
                       gene_col = "target", split_col = "CellComp",
                       file_out = "cosmx_brain/tx_spots", z = "all",
                       BPPARAM = SerialParam(progressbar = TRUE))
rowGeometries(sfe) <- spots_sf

# Test data-------------
# Just take the first FOV and 2 of the z-planes. May need to downsample spots.
spots <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_tx_file.csv")
meta <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_metadata_file.csv")
polys <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter-polygons.csv")
# I think 0 means not in any cell
mat <- fread("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_exprMat_file.csv")

# Take a small bbox from FOV 1
polys <- polys |> filter(fov == 1L)
polys_sf <- df2sf(polys, spatialCoordsNames = c("x_global_px", "y_global_px"),
                  geometryType = "POLYGON", id_col = "cellID")
all(st_is_valid(polys_sf))
bbox_use <- st_bbox(polys_sf) # 4245 px in x, 4255 in y
bbox_new <- c(xmin = bbox_use["xmin"], xmax = bbox_use["xmin"] + 425*2, 
              ymin = bbox_use["ymin"], ymax = bbox_use["ymin"] + 425*2)
names(bbox_new) <- c("xmin", "xmax", "ymin", "ymax")
bbox_sf <- st_as_sfc(st_bbox(bbox_new))

polys_sf_small <- polys_sf |> 
    st_filter(bbox_sf, .predicate = st_covered_by)
# 27 cells, fine
cells_use <- as.integer(polys_sf_small$ID)

polys <- polys |> filter(cellID %in% cells_use)
meta <- meta |> filter(fov == 1L, cell_ID %in% cells_use)
mat <- mat |> filter(fov == 1L, cell_ID %in% c(0L, cells_use))
spots <- spots |> filter(fov == 1L, z %in% 0:1, cell_ID %in% c(0L, cells_use))
# spots file still a bit too large

dir.create("inst/extdata/cosmx")
write_csv(mat, "inst/extdata/cosmx/Run5642_S3_Quarter_exprMat_file.csv")
write_csv(spots, "inst/extdata/cosmx/Run5642_S3_Quarter_tx_file.csv")
write_csv(polys, "inst/extdata/cosmx/Run5642_S3_Quarter-polygons.csv")
write_csv(meta, "inst/extdata/cosmx/Run5642_S3_Quarter_metadata_file.csv")
