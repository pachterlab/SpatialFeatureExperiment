library(vroom)
library(rhdf5)
library(EBImage)
# Here I'll use the first FOV for testing and examples
h52poly_fov <- function(fn, i) {
    l <- rhdf5::h5dump(fn)[[1]]
    cell_ids <- names(l)
    geometries <- lapply(l, function(m) sf::st_polygon(list(t(m[["zIndex_0"]]$p_0$coordinates[,,1]))))
    df <- data.frame(geometry = sf::st_sfc(geometries),
                     ID = cell_ids,
                     fov = i)
    sf::st_sf(df)
}
setwd("~/SFEData")
fns <- list.files("cell_boundaries", "*.hdf5", full.names = TRUE)
#polys <- bpmapply(h52poly_fov, fn = fns, i = seq_along(fns), SIMPLIFY = FALSE,
#                  BPPARAM = SnowParam(20, progressbar = TRUE))
#polys <- do.call(bind_rows, polys)
# Only keep the first z plane since they're all the same anyway

mat <- vroom("Liver1Slice1_cell_by_gene.csv", col_types = cols(...1 = "c")) |> as.data.frame()
mat <- mat[mat$...1 %in% polys$ID,]
polys <- polys[match(mat$...1, polys$ID),]
metadata <- vroom("Liver1Slice1_cell_metadata.csv", col_types = cols(...1 = "c")) |> as.data.frame()
metadata <- metadata[match(mat$...1, metadata$...1),]
rownames(metadata) <- metadata$...1
metadata$...1 <- NULL
rownames(mat) <- mat$...1
mat$...1 <- NULL

dir_use <- "~/SpatialFeatureExperiment/inst/extdata/vizgen"
dir.create(dir_use)
write.csv(mat, file.path(dir_use, "Liver1Slice1_cell_by_gene.csv"),
          quote = FALSE, row.names = TRUE)
write.csv(metadata, file.path(dir_use, "Liver1Slice1_cell_metadata.csv"),
          quote = FALSE, row.names = TRUE)

dir.create(file.path(dir_use, "cell_boundaries"))
ids <- rownames(mat)

new_fn <- file.path(dir_use, "cell_boundaries", "feature_data_0a.hdf5")
h5createFile(new_fn)
h5createGroup(new_fn, "featuredata")
for (i in ids) {
    gn <- file.path("featuredata", i, "zIndex_0/p_0/coordinates")
    h5createGroup(new_fn, file.path("featuredata", i))
    h5createGroup(new_fn, file.path("featuredata", i, "zIndex_0"))
    h5createGroup(new_fn, file.path("featuredata", i, "zIndex_0", "p_0"))
    coords <- h5read(fns[1], gn)
    h5write(coords, new_fn, gn)
}

l <- h5ls(fns[1])

polys <- h52poly_fov(fns[1], 1)

#file.copy("cell_boundaries/feature_data_0.hdf5",
#          file.path(dir_use, "cell_boundaries/feature_data_0.hdf5"))

# Get image
dir.create(file.path(dir_use, "images"))

library(terra)
library(tidyterra)
library(sf)
polyt <- rast("mosaic_PolyT_z0.tif")
plot(polyt)
library(jsonlite)
manifest <- read_json("manifest.json", simplifyVector = TRUE)
ext_use <- setNames(manifest$bbox_microns, c("xmin", "ymin", "xmax", "ymax"))
ext(polyt) <- ext_use[c("xmin", "xmax", "ymin", "ymax")]

bbox_orig <- bbox <- st_bbox(polys)
dist1 <- bbox["ymin"] - ext_use["ymin"]
dist2 <- bbox["ymax"] - ext_use["ymin"]
bbox["ymin"] <- ext_use["ymax"] - dist2
bbox["ymax"] <- ext_use["ymax"] - dist1

polyt_sub <- terra::crop(polyt, ext(bbox[c("xmin", "xmax", "ymin", "ymax")]))
polyt_sub <- terra::shift(polyt_sub, dy = - (ext_use["ymax"] - ext_use["ymin"] - dist2 + dist1))
polyt_sub <- terra::flip(polyt_sub)
polyt_sub <- terra::shift(polyt_sub, dy = bbox_orig["ymin"] - as.vector(ext(polyt_sub))["ymin"])

plot(polyt_sub)

st_crs(polys) <- 2154
crs(polyt_sub) <- "epsg:2154"
library(ggplot2)
# Check if they line up
ggplot() +
    geom_spatraster(data = polyt_sub) +
    geom_sf(data = polys, fill = NA, color = "white") +
    scale_fill_viridis_c() +
    coord_sf(datum = 2154)

# It worked. I'm saving this script and the image
polyt_sub2 <- terra::aggregate(polyt_sub, fact = 4)
img_out <- as.array(flip(polyt_sub2))
img_out[is.nan(img_out)] <- 0

library(tiff)
writeTIFF(img_out/max(img_out), file.path(dir_use, "images/mosaic_PolyT_z0.tif"))

# Write toy manifest.json
mnfst <- manifest
mnfst$mosaic_width_pixels <- ncol(polyt_sub2)
mnfst$mosaic_height_pixels <- nrow(polyt_sub2)
mnfst$bbox_microns <- as.vector(ext(polyt_sub2))[c("xmin", "ymin", "xmax", "ymax")] |>
    unname()
mnfst$hor_num_tiles_box <- 1
mnfst$vert_num_tiles_box <- 1
mnfst$mosaic_files <- NULL
mnfst$mosaic_pyramid_files <- NULL
write_json(mnfst, file.path(dir_use, "images/manifest.json"), pretty = TRUE,
           auto_unbox = TRUE)

# Look at the Cellpose parquet file, from a different dataset
library(sfarrow)
df <- st_read_parquet("cellpose_micron_space.parquet")
unique(df$Type)
unique(df$ZLevel)
all(is.na(df$ParentType))

# ID column has one entry for each segmentation in each z-plane
# EntityID column is for cell IDs across all z-planes

df[df$ZIndex == 3L,] |>
    ggplot() + geom_sf()
df2 <- df[df$ZIndex == 3L,] # Looks like all z-planes are the same anyway
is_multi <- st_geometry_type(df2$Geometry) == "MULTIPOLYGON"
l <- lengths(df2$Geometry)
is_multi <- which(l > 1L)

areas <- lapply(df2$Geometry[is_multi], function(x) {
    st_area(st_cast(st_sfc(x), "POLYGON"))
})
l2 <- lengths(areas)
table(l2)
summary(unlist(areas))
# Size of the larger piece
areas_large <- vapply(areas, function(x) x[which.max(x)], FUN.VALUE = numeric(1))
summary(areas_large)

# Size of the second largest piece
areas_small <- vapply(areas, function(x) sort(x, decreasing = TRUE)[2],
                      FUN.VALUE = numeric(1))
summary(areas_small)

areas_prop <- areas_small/areas_large
summary(areas_prop)
# Based on these, I think threshold of 15 micron^2 may be reasonable.
# I'll remove any cell all of whose pieces are too small.

# Back to toy example
parq <- polys
# Not using that field anyway so can be random
parq$ID <- sample(seq_len(nrow(parq)) - 1, nrow(parq))
parq$EntityID <- polys$ID
parq$ZIndex <- 0
parq$Type <- "cell"
parq$ZLevel <- 1.5
parq$ParentID <- parq$ParentType <- parq$Name <- NA
parq$X__index_level_0__ <- parq$ID
names(parq)[names(parq) == "geometry"] <- "Geometry"
st_geometry(parq) <- "Geometry"
parq <- parq[,names(df)]

st_write_parquet(parq, file.path(dir_use, "cellpose_micron_space.parquet"))
