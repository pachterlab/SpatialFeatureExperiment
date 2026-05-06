# Find Visium HD spatial neighborhood graph

Visium HD spots are arranged in a square grid. This function finds
either a rook or a queen spatial neighborhood graph for the spots.
`colData` of the SFE object must have columns `array_row` and
`array_col`.

## Usage

``` r
findVisiumHDGraph(x, style = "W", queen = FALSE, zero.policy = TRUE)
```

## Arguments

- x:

  An SFE object with Visium HD data with one sample with the required
  information in its `colData`.

- style:

  `style` can take values “W”, “B”, “C”, “U”, “minmax” and “S”

- queen:

  Logical. Default is `FALSE`, using rook neighbors.

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

## Value

A `listw` object for the graph.
