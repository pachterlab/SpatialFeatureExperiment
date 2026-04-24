# Change sample IDs

Change sample IDs in all fields of the SFE object where sample IDs are
present, not just the colData.

## Usage

``` r
changeSampleIDs(sfe, replacement)
```

## Arguments

- sfe:

  A `SpatialFeatureExperiment` object.

- replacement:

  A named character vector whose names are the existing sample IDs to be
  changed and whose values are the corresponding replacements.

## Value

An SFE object.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe <- changeSampleIDs(sfe, c(Vis5A = "sample01"))
sampleIDs(sfe)
#> [1] "sample01"
```
