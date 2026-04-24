# Get unit of a SpatialFeatureExperiment

Length units can be microns or pixels in full resolution image in SFE
objects.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
unit(x)
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

## Value

A string for the name of the unit. At present it's merely a string and
`udunits` is not used.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
SpatialFeatureExperiment::unit(sfe)
#> [1] "full_res_image_pixels"
```
