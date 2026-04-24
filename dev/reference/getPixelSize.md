# Get physical size of pixels

This function gets physical size of pixels in each resolution of a
OME-TIFF pyramid in
[`BioFormatsImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/BioFormatsImage.md).

## Usage

``` r
getPixelSize(file, resolution = 1L)
```

## Arguments

- file:

  Path to an OME-TIFF file.

- resolution:

  Which resolution to query; 1 means the highest resolution. The pixels
  will be larger for the lower resolutions.

## Value

Numeric vector of length 2 of pixel size in x and y. Usually they're the
same.

## Examples

``` r
library(SFEData)
fp <- tempfile()
dir_use <- XeniumOutput("v1", file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> The downloaded files are in /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce812969c55/xenium_lr 
# RBioFormats null pointer error
try(getPixelSize(file.path(dir_use, "morphology_focus.ome.tif")))
#> [1] 1.7 1.7
getPixelSize(file.path(dir_use, "morphology_focus.ome.tif"))
#> [1] 1.7 1.7
unlink(dir_use, recursive = TRUE)
```
