# Whether an SFE object contains out of memory data

Out of memory data, such as `DelayedArray`, some `SpatRasterImage`
objects, and `BioFormatsImage`, will break if saved as RDS. This method
of
[`containsOutOfMemoryData`](https://rdrr.io/pkg/BiocGenerics/man/containsOutOfMemoryData.html)
checks if an SFE object has out of memory data, specifically the images.
Having out of memory data will result into an error when `saveRDS` is
called; we recommend using the `alabaster.sfe` package instead.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
containsOutOfMemoryData(object)
```

## Arguments

- object:

  An SFE object

## Value

TRUE or FALSE

## Examples

``` r
outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
samples <- file.path(outdir, paste0("sample0", 1:2))
sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample01
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample02
containsOutOfMemoryData(sfe)
#> [1] TRUE
```
