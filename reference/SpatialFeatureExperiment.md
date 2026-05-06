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
  annotGeometryType = "POLYGON",
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
  of one of the following:

  An `sf` data frame

  :   The geometry column specifies the geometry of the entities.

  An ordinary data frame specifying centroids

  :   Column names for the coordinates are specified in the
      `spatialCoordsNames` argument. For Visium and ST, in addition to
      the centroid coordinate data frame, the spot diameter in the same
      unit as the coordinates can be specified in the `spotDiamter`
      argument.

  An ordinary data frame specifying polygons

  :   Also use `spatialCoordsNames`. There should an additional column
      "ID" to specify which vertices belong to which polygon. The
      coordinates should not be in list columns. Rather, the data frame
      should look like it is passed to
      [`ggplot2::geom_polygon`](https://ggplot2.tidyverse.org/reference/geom_polygon.html).
      If there are holes, then there must also be a column "subID" that
      differentiates between the outer polygon and the holes.

  In all cases, the data frame should specify the same number of
  geometries as the number of columns in the gene count matrix. If the
  column "barcode" is present, then it will be matched to column names
  of the gene count matrix. Otherwise, the geometries are assumed to be
  in the same order as columns in the gene count matrix. If the
  geometries are specified in an ordinary data frame, then it will be
  converted into `sf` internally. Named list of data frames because each
  entity can have multiple geometries, such as whole cell and nuclei
  segmentations. The geometries are assumed to be POINTs for centroids
  and POLYGONs for segmentations. If polygons are specified in an
  ordinary data frame, then anything with fewer than 3 vertices will be
  removed. For anything other than POINTs, attributes of the geometry
  will be ignored.

- rowGeometries:

  Geometry associated with genes or features, which correspond to rows
  of the gene count matrix.

- annotGeometries:

  Geometry of entities that do not correspond to columns or rows of the
  gene count matrix, such as tissue boundary and pathologist annotations
  of histological regions, and nuclei segmentation in a Visium dataset.
  Also a named list as in `colGeometries`. The ordinary data frame may
  specify POINTs, POLYGONs, or LINESTRINGs, or their MULTI versions.
  Each data frame can only specify one type of geometry. For MULTI
  versions, there must be a column "group" to identify each MULTI
  geometry.

- spotDiameter:

  Spot diameter for technologies with arrays of spots of fixed diameter
  per slide, such as Visium, ST, DBiT-seq, and slide-seq. The diameter
  must be in the same unit as the coordinates in the \*Geometry
  arguments. Ignored for geometries that are not POINT or MULTIPOINT.

- annotGeometryType:

  Character vector specifying geometry type of each element of the list
  if `annotGeometry` is specified. Each element of the vector must be
  one of POINT, LINESTRING, POLYGON, MULTIPOINT, MULTILINESTRING, and
  MULTIPOLYGON. Must be either length 1 (same for all elements of the
  list) or the same length as the list. Ignored if the corresponding
  element is an `sf` object.

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
```
