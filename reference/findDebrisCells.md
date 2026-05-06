# Identify cells in small bits outside the main piece of tissue

This is used for quality control (QC), as the small bits are likely to
be low quality technical artifacts and are not informative to spatial
analyses. Please confirm the quality of those cells by checking the
histology image with ImageJ or QuPath.

## Usage

``` r
# S4 method for class 'matrix'
findDebrisCells(
  x,
  max_cells = 5,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'sfc'
findDebrisCells(
  x,
  max_cells = 5,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'sf'
findDebrisCells(
  x,
  max_cells = 5,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)

# S4 method for class 'SpatialExperiment'
findDebrisCells(
  x,
  max_cells = 5,
  distance_cutoff = 50,
  BNPARAM = NULL,
  BPPARAM = SerialParam()
)
```

## Arguments

- x:

  Object with spatial coordinates of cells. Can be a `SpatialExperiment`
  or `SpatialFeatureExperiment` object, a matrix with 2 columns for x
  and y coordinates of cells, or a `sf` or `sfc` object with cell
  geometries.

- max_cells:

  Maximum number of cells for a clump of cells to be considered debris.

- distance_cutoff:

  Minimum distance of cell to the tissue for it to be considered debris,
  in the same unit as in `x`.

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

Depends on the method:

- Spatial(Feature)Experiment:

  The same object with a logical column "is_debris" added to `colData`.

- Matrix and `sf(c)`:

  A logical vector indicating whether each cell is debris.

## Details

How this function works: A distance-based spatial neighborhood graph is
computed, with `distance_cutoff` as the distance cutoff. Then disjoint
connected subgraphs are found. Cells in subgraphs with `max_cells` or
fewer cells are considered debris.
