# Dimension geometry methods

"Dimension geometry" refers to Simple Feature (`sf`) geometries
associated with rows (features, genes) or columns (cells or spots) of
the gene count matrix in the `SpatialFeatureExperiment` object. For each
dimension, the number of rows in the `sf` data frame specifying the
geometries must match the size of the dimension of interest. For
example, there must be the same number of rows in the `sf` data frame
describing cells as there are cells in the gene count matrix. This page
documents getters and setters for the dimension geometries. The getters
and setters are implemented in a way similar to those of `reducedDims`
in `SingleCellExperiment`.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
dimGeometries(x, MARGIN = 2, withDimnames = TRUE)

# S4 method for class 'SpatialFeatureExperiment'
dimGeometries(x, MARGIN, withDimnames = TRUE, translate = TRUE, ...) <- value

# S4 method for class 'SpatialFeatureExperiment'
dimGeometryNames(x, MARGIN)

# S4 method for class 'SpatialFeatureExperiment,numeric,character'
dimGeometryNames(x, MARGIN) <- value

# S4 method for class 'SpatialFeatureExperiment'
dimGeometry(x, type = 1L, MARGIN, sample_id = 1L, withDimnames = TRUE)

# S4 method for class 'SpatialFeatureExperiment'
dimGeometry(
  x,
  type = 1L,
  MARGIN,
  sample_id = 1L,
  withDimnames = TRUE,
  translate = TRUE,
  ...
) <- value
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- MARGIN:

  As in [`apply`](https://rdrr.io/r/base/apply.html). 1 stands for rows
  and 2 stands for columns.

- withDimnames:

  Logical. If `TRUE`, then the dimnames (colnames or rownames) of the
  gene count matrix should correspond to row names of the `sf` data
  frames of interest.

- translate:

  Logical. Only used if
  [`removeEmptySpace`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/removeEmptySpace.md)
  has been run of the SFE object. If that's the case, this argument
  indicates whether the new value to be assigned to the geometry is in
  the coordinates prior to removal of empty space so it should be
  translated to match the new coordinates after removing empty space.
  Default to `TRUE`.

- ...:

  `spatialCoordsNames, spotDiameter, geometryType` passed to
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md).
  Defaults are the same as in
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md).
  For `dimGeometries<-` only: `geometryType` can be a character vector
  of the geometry type of each data frame in the list of the same length
  as the list if the data frames specify different types of geometries.

- value:

  Value to set. For `dimGeometry`, must be a `sf` data frame with the
  same number of rows as size in the dimension of interest, or an
  ordinary data frame that can be converted to such a `sf` data frame
  (see
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md)).
  For `dimGeometries`, must be a list of such `sf` or ordinary data
  frames.

- type:

  An integer specifying the index or string specifying the name of the
  \*Geometry to query or replace. If missing, then the first item in the
  \*Geometries will be returned or replaced.

- sample_id:

  Sample ID to get or set geometries.

## Value

Getters for multiple geometries return a named list. Getters for names
return a character vector of the names. Getters for single geometries
return an `sf` data frame. Setters return an SFE object.

## See also

\[colGeometries()\], \[rowGeometries()\]

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache

# Get all column geometries as a named list
# Use MARGIN = 1 or rowGeometry/ies for rowGeometries
cgs <- dimGeometries(sfe, MARGIN = 2)
# Or equivalently
cgs <- colGeometries(sfe)

# Set all column geometries with a named list
dimGeometries(sfe, MARGIN = 2) <- cgs
# Or equivalently
colGeometries(sfe) <- cgs

# Get names of column geometries
cgns <- dimGeometryNames(sfe, MARGIN = 2)
cgns <- colGeometryNames(sfe)

# Set column geometry names
dimGeometryNames(sfe, MARGIN = 2) <- cgns
colGeometryNames(sfe) <- cgns

# Get a specific column geometry by name
spots <- dimGeometry(sfe, "spotPoly", MARGIN = 2)
spots <- colGeometry(sfe, "spotPoly")
# Or equivalently, the wrapper specifically for Visium spot polygons,
# for the name "spotPoly"
spots <- spotPoly(sfe)
# Other colGeometry wrappers for specific names:
# ROIPoly (for LCM and GeoMX DSP), cellSeg and nucSeg (for MERFISH; would
# query annotGeometries for Visium)
# rowGeometry wrappers for specific names: txSpots (MERFISH transcript spots)
# By index
spots <- colGeometry(sfe, 1L)

# Multiple samples, only get geometries for one sample
sfe2 <- McKellarMuscleData("small2")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe_combined <- cbind(sfe, sfe2)
spots1 <- colGeometry(sfe, "spotPoly", sample_id = "Vis5A")
spots2 <- spotPoly(sfe_combined, sample_id = "sample02")
# Get geometries for multiple samples
spots3 <- spotPoly(sfe_combined, sample_id = c("Vis5A", "sample02"))
# All samples
spots3 <- spotPoly(sfe_combined, sample_id = "all")

# Set specific column geometry by name
colGeometry(sfe, "foobar") <- spots
# Or use wrapper
spotPoly(sfe) <- spots
# Specify sample_id
colGeometry(sfe_combined, "foobar", sample_id = "Vis5A") <- spots1
# Only entries for the specified sample are set.
foobar <- colGeometry(sfe_combined, "foobar", sample_id = "sample02")
```
