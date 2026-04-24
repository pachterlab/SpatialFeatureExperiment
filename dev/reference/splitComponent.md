# Split by graph component

Split either a SFE object or `sf` or `sfc` by graph component of a
spatial neighborhood graph.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
splitComponent(x, colGraphName = 1L, min_cells = 100)

# S4 method for class 'sf'
splitComponent(
  x,
  min_cells = 100,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'sfc'
splitComponent(
  x,
  min_cells = 100,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)
```

## Arguments

- x:

  A SFE object, or `sf` or `sfc`.

- colGraphName:

  Name of graph to use for `splitComponent`. We recommend distance based
  neighbors (\`dnearneigh\` in
  [`findSpatialNeighbors`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/findSpatialNeighbors.md)),
  and recommend NOT using k nearest neighbors (\`knearneigh\`) or
  triangulation (\`tri2nb\`).

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

A list of the same object type as the input, each element for one
component.
