# Toy examples for df2sf
library(sf)
library(tidyverse)
# For all, check: Correct number of rows, correct geometry type, that the
# original x and y columns were removed and replaced by geometry column

# Plain points, no spotDiameter
set.seed(29)
pts <- matrix(runif(10), ncol = 2)
colnames(pts) <- c("x", "y")
pts_geometry <- st_sfc(apply(pts, 1, st_point, simplify = FALSE))
pts_df <- as.data.frame(pts)
pts_df$ID <- sample(LETTERS, 5)
pts_sf <- st_sf(ID = pts_df$ID, geometry = pts_geometry,
                sf_column_name = "geometry")
saveRDS(pts_df, "inst/testdata/pts_df.rds")
saveRDS(pts_sf, "inst/testdata/pts_sf.rds")

# Points with spotDiameter
pts_sf_dia <- st_buffer(pts_sf, dist = 0.05)
saveRDS(pts_sf_dia, "inst/testdata/pts_sf_dia.rds")

# Below: Check presence of columns ID, sample_id, and group and subID if applicable
# MULTIPOINT
multipts_df <- pts_df
multipts_df$group <- c("A", "A", "B", "B", "C")
multipts_df$sample_id <- c(rep("sample01", 2), rep("sample02", 3))
multipts_sf <- multipts_df %>%
  group_nest(group) %>%
  mutate(geometry = map(data, function(x) {
    st_multipoint(as.matrix(x[,c("x", "y")]))
  }),
  geometry = st_sfc(geometry)) %>%
  select(-data) %>%
  st_sf(sf_column_name = "geometry")
multipts_sf$sample_id <- c("sample01", "sample02", "sample02")
saveRDS(multipts_df, "inst/testdata/multipts_df.rds")
saveRDS(multipts_sf, "inst/testdata/multipts_sf.rds")

# Also expect error when the same group has more than one unique sample_ids
multipts_df_wrong_sample <- multipts_df %>%
  mutate(sample_id = c(rep("sample01", 3), rep("sample02", 2)))
saveRDS(multipts_df_wrong_sample, "inst/testdata/multipts_df_wrong_sample.rds")

# Below: Also check removal of items with too few vertices
# LINESTRING
ls_df <- pts_df %>%
  mutate(ID = c(rep("A", 2), rep("B", 3)))
ls_sf <- ls_df %>%
  group_nest(ID) %>%
  mutate(geometry = map(data, function(x) {
    st_linestring(as.matrix(x[,c("x", "y")]))
  }),
  geometry = st_sfc(geometry)) %>%
  select(-data) %>%
  st_sf(sf_column_name = "geometry")
saveRDS(ls_df, "inst/testdata/ls_df.rds")
saveRDS(ls_sf, "inst/testdata/ls_sf.rds")

# Too few vertices, expect warning message from .df2sf_check
ls_df_singleton <- ls_df %>%
  mutate(ID = c(rep("A", 2), rep("B", 2), "C"))
ls_sf_singleton <- ls_df_singleton %>%
  filter(ID != "C") %>%
  group_nest(ID) %>%
  mutate(geometry = map(data, function(x) {
    st_linestring(as.matrix(x[,c("x", "y")]))
  }),
  geometry = st_sfc(geometry)) %>%
  select(-data) %>%
  st_sf(sf_column_name = "geometry")
saveRDS(ls_df_singleton, "inst/testdata/ls_df_singleton.rds")
saveRDS(ls_sf_singleton, "inst/testdata/ls_sf_singleton.rds")

# MULTILINESTRING
multils_df <- ls_df %>%
  mutate(group = "G")

# POLYGON
# Without holes
# With holes
# Check that the geometry is POINT when the df is de facto specifying point
# regardless of the geometryType argument.

# MULTIPOLYGON
