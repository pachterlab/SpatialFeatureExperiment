# Get all unique sample IDs

The title is self-explanatory.

## Usage

``` r
sampleIDs(sfe)
```

## Arguments

- sfe:

  A `SpatialFeatureExperiment` object.

## Value

A character vector of all unique entries of the `sample_id` column in
`colData(x)`.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sampleIDs(sfe)
#> [1] "Vis5A"
```
