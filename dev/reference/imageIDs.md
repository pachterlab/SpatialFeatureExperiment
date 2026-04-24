# Show all image_ids in the SFE object

The title is self-explanatory. Some functions require `image_id` to get
or set images.

## Usage

``` r
imageIDs(sfe)
```

## Arguments

- sfe:

  A `SpatialFeatureExperiment` object.

## Value

A character vector of `image_ids`.

## Examples

``` r
fp <- system.file(file.path("extdata", "sample01"),
package = "SpatialFeatureExperiment")
sfe <- read10xVisiumSFE(fp, type = "sparse")
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample01
imageIDs(sfe)
#> [1] "lowres" "hires" 
```
