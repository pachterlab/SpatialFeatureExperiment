# Subsetting SpatialFeatureExperiment objects

The SFE method has special treatment for the spatial graphs. In `listw`,
the neighbors are indicated by indices, which will change after
subsetting. The `SFE_graph_subset` option determines whether the graphs
are subsetted or reconstructed. In the default
(`options(SFE_graph_subset = TRUE)`), the graphs are subsetted, in which
case singletons may be produced. For
`options(SFE_graph_subset = FALSE)`, which is the behavior of versions
earlier than Bioc 3.20, the graphs are reconstructed with the parameters
recorded in an attribute of the graphs. This option can result into
different graphs. For example, suppose we start with a k nearest
neighbor graph. After subsetting, cells at the boundary of the region
used to subset the SFE object may lose some of their neighbors. In
contrast, when the graph is reconstructed, these same edge cells will
gain other cells that remain after subsetting as neighbors in the new
KNN graph.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- i:

  Row indices for subsetting.

- j:

  column indices for subsetting.

- ...:

  Passed to the `SingleCellExperiment` method of `[`.

- drop:

  Only used if graphs are reconstructed
  (`options(SFE_graph_subset = FALSE)`). If `TRUE` then `colGraphs` are
  dropped but `annotGraphs` are kept.

## Value

A subsetted `SpatialFeatureExperiment` object.

## Details

The option `SFE_graph_subset` was introduced because subsetting is
usually faster than reconstructing and in some cases such as
distance-based neighbors and Visium spot adjacency give the same
results. It was introduced also because of the development of
`alabster.sfe` for a language-agnostic on-disk serialization of SFE
objects and some parameters used to construct graphs have special
classes whose `alabaster` methods have not been implemented, such as
`BPPARAM` and `BNPARAM`, so when reconstructing, the defaults for those
arguments will be used.

The edge weights will be recomputed from the binary neighborhood
indicator with the same normalization style as the original graph, such
as "W" for row normalization. When distance-based edge weights are used
instead of the binary indicator, the edge weights will be re-normalized,
which is mostly some rescaling. This should give the same results as
recomputing the distance based edge weights for styles "raw", "W", and
"B" since the distances themselves don't change, but the effects of
other more complicated styles of re-normalization on spatial statistics
should be further investigated.

By default, upon subsetting, the images are cropped to the bounding box
of the remaining cells. However, when the image is large and the
bounding box contains most of the original image, cropping is slow.
Cropping can be disabled by `options(SFE_subset_crop = FALSE)`. Also,
when the remaining part of the image is larger than a threshold, the
image will not be cropped; the threshold can be set with the
`SFE_subset_crop_max` option, such as
`options(SFE_subset_crop_max = "100MB")`.

## Examples

``` r
# Just like subsetting matrices and SingleCellExperiment
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe_subset <- sfe[seq_len(10), seq_len(10), drop = TRUE]
# Gives warning as graph reconstruction fails
# \donttest{
sfe_subset <- sfe[seq_len(10), seq_len(10)]
# }
```
