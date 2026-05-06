# Find spatial neighborhood graph

This function wraps all spatial neighborhood graphs implemented in the
package `spdep` for the `SpatialFeatureExperiment` (SFE) class, to find
spatial neighborhood graphs for the entities represented by columns or
rows of the gene count matrix in the SFE object or spatial entities in
the `annotGeometries` field of the SFE object. Results are stored as
`listw` objects in the `spatialGraphs` field of the SFE object, as
`listw` is used in many methods that facilitate the spatial neighborhood
graph in the `spdep`, `spatialreg`, and `adespatial`. The edge weights
of the graph in the `listw` object are by default style W (see
[`nb2listw`](https://r-spatial.github.io/spdep/reference/nb2listw.html))
and the unweighted neighbor list is in the `neighbours` field of the
`listw` object.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
findSpatialNeighbors(
  x,
  sample_id = "all",
  type = "spatialCoords",
  MARGIN = 2,
  method = c("tri2nb", "knearneigh", "dnearneigh", "gabrielneigh", "relativeneigh",
    "soi.graph", "poly2nb"),
  dist_type = c("none", "idw", "exp", "dpd"),
  glist = NULL,
  style = c("raw", "W", "B", "C", "U", "minmax", "S"),
  nn_method = c("bioc", "spdep"),
  alpha = 1,
  dmax = NULL,
  BPPARAM = SerialParam(),
  BNPARAM = KmknnParam(),
  zero.policy = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`SpatialFeatureExperiment`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment.md)
  object.

- sample_id:

  Which sample(s) in the SFE object to use for the graph. Can also be
  "all", which means this function will compute the graph for all
  samples independently.

- type:

  Name of the geometry associated with the MARGIN of interest for which
  to compute the graph.

- MARGIN:

  Just like in [`apply`](https://rdrr.io/r/base/apply.html), where 1
  stands for row, 2 stands for column. Here, in addition, 3 stands for
  annotation, to query the
  [`annotGeometries`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md),
  such as nuclei segmentation in a Visium data

- method:

  Name of function in the package `spdep` to use to find the spatial
  neighborhood graph.

- dist_type:

  Type of distance-based weight. "none" means not using distance-based
  weights; the edge weights of the spatial neighborhood graph will be
  entirely determined by the `style` argument. "idw" means inverse
  distance weighting. "exp" means exponential decay. "dpd" means
  double-power distance weights. See
  [`nb2listwdist`](https://r-spatial.github.io/spdep/reference/nb2listwdist.html)
  for details.

- glist:

  list of general weights corresponding to neighbours

- style:

  `style` can take values “W”, “B”, “C”, “U”, “minmax” and “S”

- nn_method:

  Method to find k nearest neighbors and distance based neighbors. Can
  be either "bioc" or "spdep". For "bioc", methods from `BiocNeighbors`
  are used. For "spdep", methods from the `spdep` package are used. The
  "bioc" option is more scalable to larger datasets and supports
  multithreading.

- alpha:

  Only relevant when `dist_type = "dpd"`.

- dmax:

  Only relevant when `dist_type = "dpd"`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object for multithreading. Only used for k nearest neighbor and
  distance based neighbor with `nn_method = "bioc"`.

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

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

- ...:

  Extra arguments passed to the `spdep` function stated in the `method`
  argument, such as `k`, `use_kd_tree`, `d1`, `d2`, `nnmult`, `sym`, and
  `quadsegs`. Note that any arguments about using longitude and
  latitude, which are irrelevant, are ignored.

## Value

For one sample, then a `listw` object representing the graph, with an
attribute "method" recording the function used to build the graph, its
arguments, and information about the geometry for which the graph was
built. The attribute is used to reconstruct the graphs when the SFE
object is subsetted since some nodes in the graph will no longer be
present. If sample_id = "all" or has length \> 1, then a named list of
`listw` objects, whose names are the sample_ids. To add the list for
multiple samples to a SFE object, specify the `name` argument in the
[`spatialGraphs`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
replacement method, so graph of the same name will be added to the SFE
object for each sample.

## Note

`style = "raw"` is only applicable when `dist_type` is not "none". If
`dist_type = "none"` and `style = "raw"`, then style will default to
"W". Using distance based weights does not supplant finding a spatial
neighborhood graph. The spatial neighborhood graph is first found and
then its edges weighted based on distance in this function.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
# sample_id is optional when only one sample is present
g <- findSpatialNeighbors(sfe, sample_id = "Vis5A")
attr(g, "method")
#> $FUN
#> [1] "findSpatialNeighbors"
#> 
#> $package
#> $package[[1]]
#> [1] "SpatialFeatureExperiment"
#> 
#> $package[[2]]
#> [1] "1.14.0"
#> 
#> 
#> $args
#> $args$method
#> [1] "tri2nb"
#> 
#> $args$dist_type
#> [1] "none"
#> 
#> $args$style
#> [1] "W"
#> 
#> $args$alpha
#> [1] 1
#> 
#> $args$zero.policy
#> [1] TRUE
#> 
#> $args$sample_id
#> [1] "Vis5A"
#> 
#> $args$type
#> [1] "spatialCoords"
#> 
#> $args$MARGIN
#> [1] 2
#> 
#> 
# Returns named list for multiple samples
sfe2 <- McKellarMuscleData(dataset = "small2")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe_combined <- cbind(sfe, sfe2)
gs <- findSpatialNeighbors(sfe, sample_id = "all")
```
