# From ordinary data frame to sf to construct SFE object

While the `SpatialFeatureExperiment` constructor and `*Geometry`
replacement methods can convert properly formatted ordinary data frames
into `sf` objects which are used to store the geometries internally, the
user might want to do the conversion, check if the geometry is valid,
and inspect and fix any invalid geometries.

## Usage

``` r
df2sf(
  df,
  spatialCoordsNames = c("x", "y"),
  spotDiameter = NA,
  geometryType = c("POINT", "LINESTRING", "POLYGON", "MULTIPOINT", "MULTILINESTRING",
    "MULTIPOLYGON"),
  group_col = "group",
  id_col = "ID",
  subid_col = "subID",
  check = TRUE,
  ...
)
```

## Arguments

- df:

  An ordinary data frame, i.e. not `sf`. Or a matrix that can be
  converted to a data frame.

- spatialCoordsNames:

  Column names in `df` that specify spatial coordinates.

- spotDiameter:

  Spot diameter for technologies with arrays of spots of fixed diameter
  per slide, such as Visium, ST, DBiT-seq, and slide-seq. The diameter
  must be in the same unit as the coordinates in the \*Geometry
  arguments. Ignored for geometries that are not POINT or MULTIPOINT.

- geometryType:

  Type of geometry to convert the ordinary data frame to. If the
  geometry in `df` is de facto points, then this argument will be
  ignored and the returned `sf` will have geometry type POINT.

- group_col:

  Column to indicate which coordinates for which MULTI geometry, such as
  to identify which MULTIPOLYGON or MULTIPOINT.

- id_col:

  Column to indicate coordinates for which geometry, within a MULTI
  geometry if applicable, such as to identify which POLYGON or which
  polygon within a MULTIPOLYGON.

- subid_col:

  Column to indicate coordinates for holes in polygons.

- check:

  Logical, whether to check the input data frame for issues related to
  constructing the geometry of interese such as number of vertices per
  geometry. If `FALSE`, it will save a bit of time, which is useful when
  the input is already known to be good.

- ...:

  Other arguments passed to \`sf::st_buffer\`, mainly to make polygon
  shapes, eg Visium spot \`endCapStyle = "ROUND"\` and VisiumHD bin
  \`endCapStyle = "SQUARE"\`

## Value

An `sf` object.

## Examples

``` r
# Points, use spotDiameter to convert to circle polygons
# This is done to Visium spots
pts_df <- readRDS(system.file("extdata/pts_df.rds",
    package = "SpatialFeatureExperiment"
))
sf_use <- df2sf(pts_df, geometryType = "POINT", spotDiameter = 0.1)
# Linestring
ls_df <- readRDS(system.file("extdata/ls_df.rds",
    package = "SpatialFeatureExperiment"
))
sf_use <- df2sf(ls_df, geometryType = "LINESTRING")
# Polygon
pol_df <- readRDS(system.file("extdata/pol_df.rds",
    package = "SpatialFeatureExperiment"
))
sf_use <- df2sf(pol_df,
    geometryType = "POLYGON",
    spatialCoordsNames = c("V1", "V2")
)
# Multipolygon
mpol_df <- readRDS(system.file("extdata/mpol_df.rds",
    package = "SpatialFeatureExperiment"
))
sf_use <- df2sf(mpol_df,
    geometryType = "MULTIPOLYGON",
    spatialCoordsNames = c("V1", "V2")
)
# Multiple sample_ids present
multipts_df <- readRDS(system.file("extdata/multipts_df.rds",
    package = "SpatialFeatureExperiment"
))
sf_use <- df2sf(multipts_df, geometryType = "MULTIPOINT")
```
