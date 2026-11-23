# Constructor of SpatialFeatureExperiment object

Create a `SpatialFeatureExperiment` object.

## Usage

``` r
SpatialFeatureExperiment(
  assays,
  colData = DataFrame(),
  rowData = NULL,
  sample_id = "sample01",
  spatialCoordsNames = c("x", "y"),
  spatialCoords = NULL,
  colGeometries = NULL,
  rowGeometries = NULL,
  annotGeometries = NULL,
  spotDiameter = NA_real_,
  annotGeometryType = deprecated(),
  spatialGraphs = NULL,
  unit = c("full_res_image_pixel", "micron"),
  ...
)
```

## Arguments

- assays:

  A `list` or `SimpleList` of matrix-like elements, or a matrix-like
  object (e.g. an ordinary matrix, a data frame, a
  [DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  object from the S4Vectors package, a
  [SparseMatrix](https://rdrr.io/pkg/SparseArray/man/SparseArray-class.html)
  derivative from the SparseArray package, a
  [sparseMatrix](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  derivative from the Matrix package, a
  [DelayedMatrix](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  object from the DelayedArray package, etc...). All elements of the
  list must have the same dimensions, and dimension names (if present)
  must be consistent across elements and with the row names of
  `rowRanges` and `colData`.

- colData:

  An optional DataFrame describing the samples. Row names on `colData`,
  if present, become the column names of the returned object.

- rowData:

  `NULL` (the default) or a
  [DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  object describing the rows. Row names, if present, become the row
  names of the constructed SummarizedExperiment object. The number of
  rows of the
  [DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  must equal the number of rows of the matrices in `assays`.

- sample_id:

  A `character` sample identifier, which matches the `sample_id` in
  [`imgData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html).
  The `sample_id` will also be stored in a new column in
  [`colData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-colData.html),
  if not already present. Default = `sample01`.

- spatialCoordsNames:

  A `character` vector of column names if `*Geometries` arguments have
  ordinary data frames, to identify the columns in the ordinary data
  frames that specify the spatial coordinates. If `colGeometries` is not
  specified, then this argument will behave as in
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html),
  but `colGeometries` will be given precedence if provided.

- spatialCoords:

  A numeric matrix containing columns of spatial coordinates, as in
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).
  The coordinates are centroids of the entities represented by the
  columns of the gene count matrix. If `colGeometries` is also
  specified, then it will be given priority and a warning is issued.
  Otherwise, the `sf` representation of the centroids will be stored in
  the `colGeometry` called `centroids`.

- colGeometries:

  Geometry of the entities that correspond to the columns of the gene
  count matrix, such as cells and Visium spots. It must be a named list
  or vector of one or more of the following:

  An `sf` data frame

  :   The geometry column specifies the geometry of the entities.

  File names of GeoParquet files

  :   These will be used to represent geometries out of memory with
      `DuckDBSpatial`.

  DuckDBDataFrame with geometries

  :   As read from GeoParquet file with
      [`readGeoParquet`](https://rdrr.io/pkg/DuckDBSpatial/man/readGeoParquet.html)

  An ordinary data frame specifying centroids or polygons

  :   Deprecated

- rowGeometries:

  Geometry associated with genes or features, which correspond to rows
  of the gene count matrix.

- annotGeometries:

  Geometry of entities that do not correspond to columns or rows of the
  gene count matrix, such as tissue boundary and pathologist annotations
  of histological regions, and nuclei segmentation in a Visium dataset.
  Also a named list as in `colGeometries`.

- spotDiameter:

  Spot diameter for technologies with arrays of spots of fixed diameter
  per slide, such as Visium, ST, DBiT-seq, and slide-seq. The diameter
  must be in the same unit as the coordinates in the \*Geometry
  arguments. Ignored for geometries that are not POINT or MULTIPOINT.

- annotGeometryType:

  Deprecated. Ignored if the corresponding element is an `sf` object or
  file path of a GeoParquet file.

- spatialGraphs:

  A named list of `listw` objects (see `spdep`) for spatial neighborhood
  graphs.

- unit:

  Unit the coordinates are in, either microns or pixels in full
  resolution image.

- ...:

  Additional arguments passed to the
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  and
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  constructors.

## Value

A SFE object. If neither `colGeometries` nor `spotDiameter` is
specified, then a `colGeometry` called "centroids" will be made, which
is essentially the spatial coordinates as sf POINTs. If `spotDiameter`
is specified, but not `colGeometries`, then the spatial coordinates will
be buffered by half the diameter to get spots with the desired diameter,
and the resulting `colGeometry` will be called "spotPoly", for which
there's a convenience getter and setter,
[`spotPoly`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md).

## Examples

``` r
library(Matrix)
#> 
#> Attaching package: ‘Matrix’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     expand
data("visium_row_col")
coords1 <- visium_row_col[visium_row_col$col < 6 & visium_row_col$row < 6, ]
coords1$row <- coords1$row * sqrt(3)
cg <- df2sf(coords1[, c("col", "row")], c("col", "row"), spotDiameter = 0.7)

set.seed(29)
col_inds <- sample(seq_len(13), 13)
row_inds <- sample(seq_len(5), 13, replace = TRUE)
values <- sample(seq_len(5), 13, replace = TRUE)
mat <- sparseMatrix(i = row_inds, j = col_inds, x = values)
colnames(mat) <- coords1$barcode
rownames(mat) <- sample(LETTERS, 5)
rownames(cg) <- colnames(mat)

sfe <- SpatialFeatureExperiment(list(counts = mat),
    colData = coords1,
    spatialCoordsNames = c("col", "row"),
    spotDiameter = 0.7
)
sfe2 <- SpatialFeatureExperiment(list(counts = mat),
    colGeometries = list(foo = cg)
)

# Use DuckDBSpatial to not load geometry into memory
# Geometry is written to disk first in this toy example
library(DuckDBSpatial)
#> Loading required package: DuckDBDataFrame
#> Loading required package: bit64
#> ********************************************************
#> R-core is collecting use cases for 64-bit integers as they explore native support for these vectors.
#> 
#> See https://stat.ethz.ch/pipermail/r-devel/2026-July/084631.html and reach out to Luke Tierney (luke-tierney@uiowa.edu).
#> ********************************************************
#> 
#> Attaching package: ‘bit64’
#> The following objects are masked from ‘package:Matrix’:
#> 
#>     colSums, rowSums
#> The following objects are masked from ‘package:SummarizedExperiment’:
#> 
#>     match, order, rank
#> The following object is masked from ‘package:Biobase’:
#> 
#>     cache
#> The following objects are masked from ‘package:GenomicRanges’:
#> 
#>     as.factor, intersect, match, order, rank, setdiff, union
#> The following object is masked from ‘package:Seqinfo’:
#> 
#>     intersect
#> The following objects are masked from ‘package:IRanges’:
#> 
#>     as.factor, intersect, match, order, setdiff, table, union
#> The following objects are masked from ‘package:S4Vectors’:
#> 
#>     %in%, as.factor, intersect, match, order, rank, setdiff, setequal,
#>     union
#> The following objects are masked from ‘package:BiocGenerics’:
#> 
#>     %in%, intersect, match, order, rank, setdiff, setequal, table,
#>     union
#> The following objects are masked from ‘package:generics’:
#> 
#>     as.factor, as.ordered, intersect, is.element, setdiff, setequal,
#>     union
#> The following objects are masked from ‘package:MatrixGenerics’:
#> 
#>     colSums, rowSums
#> The following object is masked from ‘package:utils’:
#> 
#>     hashtab
#> The following objects are masked from ‘package:base’:
#> 
#>     %in%, :, array, as.factor, as.ordered, colSums, factor, intersect,
#>     is.double, is.element, match, matrix, order, rank, rowSums,
#>     setdiff, setequal, table, union
#> Loading required package: sf
#> Linking to GEOS 3.14.1, GDAL 3.12.2, PROJ 9.7.1; sf_use_s2() is TRUE
#> 
#> Attaching package: ‘DuckDBSpatial’
#> The following objects are masked from ‘package:sf’:
#> 
#>     st_contains, st_contains_properly, st_covered_by, st_covers,
#>     st_crosses, st_dimension, st_disjoint, st_distance, st_equals,
#>     st_geometry_type, st_inscribed_circle, st_is_empty, st_is_simple,
#>     st_is_within_distance, st_length, st_line_interpolate,
#>     st_line_project, st_overlaps, st_perimeter, st_touches, st_within
fn <- tempfile(fileext = ".parquet") # a file name
writeGeoParquet(cg, fn)
sfe <- SpatialFeatureExperiment(list(counts = mat),
                                colGeometries = c(foo = fn)
)
# Or use DuckDBDataFrame
cg2 <- readGeoParquet(fn)
sfe <- SpatialFeatureExperiment(list(counts = mat),
                                colGeometries = c(foo = cg2)
)

file.remove(fn)
#> [1] TRUE
```
