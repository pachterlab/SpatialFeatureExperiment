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
multils_sf_mats <- ls_df %>%
  group_nest(ID) %>%
  mutate(geometry = map(data, function(x) {
    as.matrix(x[,c("x", "y")])
  }))
multils_sf <- st_sf(ID = "G", geometry = st_sfc(st_multilinestring(multils_sf_mats$geometry)),
                    sf_column_name = "geometry")
saveRDS(multils_df, "inst/testdata/multils_df.rds")
saveRDS(multils_sf, "inst/testdata/multils_sf.rds")

# POLYGON
# Just copied this example from the first sf vignette
p1 <- rbind(c(0,0), c(1,0), c(3,2), c(2,4), c(1,4), c(0,0))
p2 <- rbind(c(1,1), c(1,2), c(2,2), c(1,1))
pol <-st_polygon(list(p1,p2))
pol_df <- as.data.frame(rbind(p1[-6,], p2[-4,]))
pol_df$ID <- "A"
pol_df$subID <- c(rep("B", 5), rep("C", 3))
pol_sf <- st_sf(ID = "A", geometry = st_sfc(pol))
saveRDS(pol_df, "inst/testdata/pol_df.rds")
saveRDS(pol_sf, "inst/testdata/pol_sf.rds")

# Check that the geometry is POINT when the df is de facto specifying point
# regardless of the geometryType argument. Don't need extra toy examples for that.

# MULTIPOLYGON
p5 <- rbind(c(3,3), c(4,2), c(4,3), c(3,3))
mpol <- st_multipolygon(list(list(p1,p2), list(p5)))
mpol_df <- as.data.frame(rbind(p1[-6,], p2[-4,], p5[-4,]))
mpol_df$ID <- c(rep("A", 8), rep("B", 3))
mpol_df$subID <- c(rep("C", 5), rep("D", 3), rep("E", 3))
mpol_df$group <- "F"
mpol_sf <- st_sf(ID = "F", geometry = st_sfc(mpol))
saveRDS(mpol_df, "inst/testdata/mpol_df.rds")
saveRDS(mpol_sf, "inst/testdata/mpol_sf.rds")
