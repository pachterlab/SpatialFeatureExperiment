# Concatenate SpatialFeatureExperiment objects

On top of the `cbind` method of `SpatialExperiment`, this method is
needed to properly merge the `spatialGraphs` field in the different SFE
objects. `rowGeometries` and `annotGeometries` also need to be combined
properly.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
cbind(..., deparse.level = 1)
```

## Arguments

- ...:

  SFE objects to cbind.

- deparse.level:

  See `?`[`rbind`](https://rdrr.io/r/base/cbind.html).

## Value

A combined SFE object.

## Examples

``` r
library(SFEData)
sfe_small <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe_small2 <- McKellarMuscleData(dataset = "small2")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
sfe2 <- cbind(sfe_small, sfe_small2)
```
