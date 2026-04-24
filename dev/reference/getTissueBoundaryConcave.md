# Get tissue boundary from concave hull of cell geometries

The concave hull will be a smoothed outline, and the smoothness can be
adjusted with the `ratio` parameter. Run
[`findDebrisCells`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/findDebrisCells.md)
to remove small bits outside the main piece tissue before running this
function.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
getTissueBoundaryConcave(
  x,
  sample_id = NULL,
  colGeometryName = 1L,
  ratio = 0.01,
  allow_holes = TRUE,
  multiple_pieces = FALSE,
  min_cells = 100,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'sfc'
getTissueBoundaryConcave(
  x,
  ratio = 0.01,
  allow_holes = TRUE,
  multiple_pieces = FALSE,
  min_cells = 100,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'sf'
getTissueBoundaryConcave(
  x,
  ratio = 0.01,
  allow_holes = TRUE,
  multiple_pieces = FALSE,
  min_cells = 100,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)
```

## Arguments

- x:

  object of class `sfg`, `sfc` or `sf`

- sample_id:

  Sample id(s) whose tissue boundaries are to be found.

- colGeometryName:

  Name of the `colGeometry` to use to infer the concave hull.

- ratio:

  numeric; fraction convex: 1 returns the convex hulls, 0 maximally
  concave hulls

- allow_holes:

  logical; if `TRUE`, the resulting concave hull may have holes

- multiple_pieces:

  Logical, whether there are multiple pieces of tissue in the same
  sample.

- min_cells:

  Minimum number of cells per graph component; components with fewer
  than this number of cells are considered debris and removed.

- distance_cutoff:

  For the `sf` and `sfc` methods, the distance cutoff used to construct
  a distance-based spatial neighborhood graph, in the same unit as the
  spatial coordinates. Anything within this distance is considered a
  neighbor.

- BNPARAM:

  A
  [`BiocNeighborParam`](https://rdrr.io/pkg/BiocNeighbors/man/BiocNeighborParam.html)
  object specifying the algorithm to find k nearest neighbors and
  distance based neighbors with `nn_method = "bioc"`. For distance based
  neighbors, only
  [`KmknnParam`](https://rdrr.io/pkg/BiocNeighbors/man/KmknnParam.html)
  and
  [`VptreeParam`](https://rdrr.io/pkg/BiocNeighbors/man/VptreeParam.html)
  are applicable.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object for multithreading. Only used for k nearest neighbor and
  distance based neighbor with `nn_method = "bioc"`.

## Value

A `sf` data frame with columns `sample_id` and `geometry`.
