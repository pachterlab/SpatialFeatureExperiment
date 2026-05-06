# Aggregate data in SFE using geometry

Gene expression and numeric columns of `colData` will be aggregated with
the function specified in `FUN`, according to another geometry supplied
and a geometry predicate (such as `st_intersects`). For example, when
the predicate is `st_intersects` and a spatial grid is used to
aggregate, then the data associated with all cells that intersect with
each grid cell will be aggregated with `FUN`, such as `mean` or `sum`.
The categorical columns will be collected into list columns, and logical
columns will be converted into numeric before applying `FUN`.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
aggregate(
  x,
  by = NULL,
  FUN = sum,
  sample_id = "all",
  colGeometryName = 1L,
  rowGeometryName = NULL,
  cellsize = NULL,
  square = TRUE,
  flat_topped = FALSE,
  new_geometry_name = "bins",
  join = st_intersects,
  sparse = FALSE,
  BPPARAM = SerialParam()
)
```

## Arguments

- x:

  An SFE object to be aggregated.

- by:

  A `sf` data frame whose geometry column is used for aggregation or
  `sfc` or for multiple samples a list of `sfc` whose names are the
  sample IDs. For multiple samples, the `sf` data frame must have a
  column `sample_id` to indicate which geometry for which sample. This
  argument is optional if `cellsize` is specified.

- FUN:

  Function to aggregate the numerical columns in `colData` and the gene
  count matrix. This can be `sum`, `mean`, or any function that takes a
  numeric matrix as input and returns a numeric vector whose length is
  same as the number of rows in the input matrix, such as `rowMedians`.
  See package `matrixStats`. Depending on the function used for
  aggregation, numeric columns of `colData` may need to be interpreted
  differently after aggregation. Aggregation is not done when
  aggregating by transcript spots in `rowGeometry`. When it's sum or
  mean, matrix multiplication is used for aggregation rather than
  calling the sum or mean function itself; this is much faster than
  looping through the bins and calling the function on each of them.

- sample_id:

  Which samples to aggregate, defaults to "all".

- colGeometryName:

  Which `colGeometry` to spatially aggregate the data, by default the
  first one.

- rowGeometryName:

  Which `rowGeometry` to spatially aggregate

- cellsize:

  numeric of length 1 or 2 with target cellsize: for square or
  rectangular cells the width and height, for hexagonal cells the
  distance between opposite edges (edge length is cellsize/sqrt(3)). A
  length units object can be passed, or an area unit object with area
  size of the square or hexagonal cell.

- square:

  logical; if `FALSE`, create hexagonal grid

- flat_topped:

  logical; if `TRUE` generate flat topped hexagons, else generate pointy
  topped

- new_geometry_name:

  Name to give to the new `colGeometry` in the output. Defaults to
  "bins".

- join:

  logical spatial predicate function to use if `by` is a simple features
  object or geometry; see
  [st_join](https://r-spatial.github.io/sf/reference/st_join.html)

- sparse:

  Logical, whether the gene count matrix from aggregating transcript
  spots should be sparse. When the bins are large, the matrix will not
  be very sparse so using sparse matrix will not save memory, but when
  the bins are small, sparsity is worth it.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object specifying parallel computing when aggregating data with
  functions other than sum and mean when aggregating cells. When
  aggregating transcript spots, this specifies parallel computing over
  genes. Defaults to `SerialParam()`.

## Value

An SFE object with `colGeometry` the same as the geometry specified in
`by` or same as the grid specified in `cellsize`. `rowGeometries` and
`rowData` remain the same as in the input `x`. `reducedDims`,
`localResults`, `colFeatureData` (and its `colGeometry`,
`annotGeometry`, and `reducedDim` counterparts), and `spatialGraphs` are
dropped because those results no longer apply after aggregation.

## Details

For smFISH-based data where the transcript spots are available, the
transcript spots can be used instead of cells to aggregate the gene
count matrix, in which case all assays other than `counts` will be
dropped and `FUN` only applies to `colData` because the transcript spots
are simply counted.

What this function does is similar to
[SEraster](https://github.com/JEFworks-Lab/SEraster) but more general
because any geometry and more aggregation function can be used, not just
regular grids, and the aggregation can be performed on the transcript
spots.

## Note

For developers: When debugging this function after calling
`devtools::load_all(".")`, you may get an error that comes from S3
dispatch of `aggregate.Vector` from the `S4Vectors` package. When that
happens, either restart the R session, or run
`setGeneric("aggregate", function(x, ...) standardGeneric("aggregate"))`
in the console to make an S4 generic as done in the `terra` package to
prioritize S4 dispatch.

## Examples

``` r
# example code
```
