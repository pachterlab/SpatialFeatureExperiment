# Get parameters used in spatial methods

The `getParams` function allows users to access the parameters used to
compute the results that may be stored in
[`colFeatureData`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/colFeatureData.md).

## Usage

``` r
getParams(
  sfe,
  name,
  local = FALSE,
  colData = FALSE,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  reducedDimName = NULL
)
```

## Arguments

- sfe:

  A `SpatialFeatureExperiment` object.

- name:

  Name used to store the results.

- local:

  Logical, whether the results of interest come from a local spatial
  method.

- colData:

  Logical, whether the results were computed for a column of
  `colData(sfe)`.

- colGeometryName:

  To get results for a `colGeometry`.

- annotGeometryName:

  To get results for an `annotGeometry`; `colGeometry` has precedence so
  this argument is ignored if `colGeometryName` is specified.

- reducedDimName:

  Name of a dimension reduction, can be seen in
  [`reducedDimNames`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).
  `colGeometryName` and `annotGeometryName` have precedence over
  `reducedDimName`.

## Value

A named list showing the parameters

## Examples

``` r
library(SFEData)
library(scater)
#> Loading required package: scuttle
#> Loading required package: ggplot2
#> 
#> Attaching package: ‘ggplot2’
#> The following object is masked from ‘package:SpatialFeatureExperiment’:
#> 
#>     unit
library(Voyager)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
colGraph(sfe, "visium") <- findVisiumGraph(sfe)
sfe <- colDataMoransI(sfe, "nCounts")
getParams(sfe, "moran", colData = TRUE)
#> $name
#> [1] "moran"
#> 
#> $package
#> [1] "spdep"
#> 
#> $version
#> [1] "1.4.2"
#> 
#> $zero.policy
#> NULL
#> 
#> $include_self
#> [1] FALSE
#> 
#> $graph_params
#> $graph_params$FUN
#> [1] "findVisiumGraph"
#> 
#> $graph_params$package
#> $graph_params$package[[1]]
#> [1] "SpatialFeatureExperiment"
#> 
#> $graph_params$package[[2]]
#> [1] "1.13.2"
#> 
#> 
#> $graph_params$args
#> $graph_params$args$style
#> [1] "W"
#> 
#> $graph_params$args$zero.policy
#> NULL
#> 
#> $graph_params$args$sample_id
#> [1] "Vis5A"
#> 
#> 
#> 
```
