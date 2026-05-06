# Update a SpatialFeatureExperiment object

Update a
[SpatialFeatureExperiment](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-class.md)
to the latest version of object structure. This is usually called by
internal functions.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
updateObject(object, ..., verbose = FALSE)

SFEVersion(object)
```

## Arguments

- object:

  An old
  [SpatialFeatureExperiment](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-class.md)
  object.

- ...:

  Additional arguments that are ignored.

- verbose:

  Logical scalar indicating whether a message should be emitted as the
  object is updated.

## Value

An updated version of `object`.

## Details

Version 1.1.4 adds package version to the SFE object. We are considering
an overhaul of the `spatialGraphs` slot in a future version using the
`sfdep` package to decouple the adjacency graph from the edge weights.

## See also

[`objectVersion`](https://rdrr.io/pkg/SingleCellExperiment/man/miscellaneous.html),
which is used to determine if the object is up-to-date.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
# First version of SFE object doesn't log SFE package version, so should be NULL
SFEVersion(sfe)
#> [1] ‘1.14.0’
sfe <- updateObject(sfe)
# See current version
SFEVersion(sfe)
#> [1] ‘1.14.0’
```
