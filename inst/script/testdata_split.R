# Make the pieces annotations to unit test split functions
# The pieces were annotated in QuPath but somehow the polygons saved from QuPath
# have a different scale from the image and geometries. Not sure if QuPath didn't
# read the physical pixel size properly or I didn't write it properly.
library(sf)
library(ggplot2)
fn <- "xenium2"
try(sfe <- readXenium(fn))
sfe <- readXenium(fn)

pieces <- st_read(system.file("extdata/pieces.geojson", package = "SpatialFeatureExperiment"),
                  drivers = "GeoJSON")
pieces <- st_geometry(pieces)
pieces <- pieces * matrix(c(1,0,0,-1), ncol = 2)
# somehow not scaled correctly though seems to be in microns
bbox_pieces <- st_bbox(pieces)
bbox_img <- ext(getImg(sfe))
sfct <- bbox_img["xmax"]/bbox_pieces["xmax"]
pieces <- pieces * sfct

# Check if it matches
ggplot() + geom_sf(data = pieces) + geom_sf(data = cellSeg(sfe)) # works

# Do the same with the other annotations
pieces2 <- st_read(system.file("extdata/subpieces.geojson", package = "SpatialFeatureExperiment"),
                   drivers = "GeoJSON", crs = NA)
pieces2$sample_id <- ifelse(pieces2$name == "region1", "sample01", "sample02")
st_geometry(pieces2) <- st_geometry(pieces2) * matrix(c(1,0,0,-1), ncol = 2)
# somehow not scaled correctly though seems to be in microns
bbox_pieces2 <- st_bbox(pieces2)
bbox_img <- ext(getImg(sfe))
sfct <- bbox_img["xmax"]/bbox_pieces2["xmax"] # should be the same for all annotations but just in case
st_geometry(pieces2) <- st_geometry(pieces2) * sfct

cont <- st_read(system.file("extdata/contiguity.geojson", package = "SpatialFeatureExperiment"),
                drivers = "GeoJSON", crs = NA)
st_geometry(cont) <- st_geometry(cont) * matrix(c(1,0,0,-1), ncol = 2)
# somehow not scaled correctly though seems to be in microns
bbox_cont <- st_bbox(cont)
bbox_img <- ext(getImg(sfe))
sfct <- bbox_img["xmax"]/bbox_cont["xmax"]
st_geometry(cont) <- st_geometry(cont) * sfct
cont$sample_id <- cont$name

# Save the annotations
saveRDS(pieces, "inst/extdata/pieces.rds")
saveRDS(pieces2[,c("sample_id", "geometry")], "inst/extdata/subpieces.rds")
saveRDS(cont[,c("sample_id", "geometry")], "inst/extdata/contiguity.rds")
