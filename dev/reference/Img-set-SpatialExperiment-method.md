# Image setter

Modify or replace images stored in a `SpatialExperiment` object. This is
different from
[`addImg`](https://rdrr.io/pkg/SpatialExperiment/man/imgData-methods.html)
which adds the image from files and can't replace existing images, which
is there to be consistent with `SpatialExperiment`. This setter here can
replace existing images with another object that inherits from
`VirtualSpatialImage`, including
[`SpatRasterImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/SpatRasterImage.md),
[`BioFormatsImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/BioFormatsImage.md),
and
[`ExtImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ExtImage.md).

## Usage

``` r
# S4 method for class 'SpatialExperiment'
Img(x, sample_id = 1L, image_id, scale_fct = 1) <- value
```

## Arguments

- x:

  A `SpatialExperiment` object, which includes SFE.

- sample_id:

  Which sample the image is associated with. Use
  [`sampleIDs`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/sampleIDs.md)
  to get sample IDs present in the SFE object.

- image_id:

  Image ID, such as "lowres" and "hires" for Visium data and "DAPI" and
  "PolyT" for Vizgen MERFISH data.

- scale_fct:

  Scale factor to convert pixels in lower resolution to those in the
  full resolution. Only relevant to image classes implemented in
  `SpatialExperiment` but not `SpatialFeatureExperiment` because the
  spatial extent of images in SFE takes precedence.

- value:

  New version of image to add, must inherit from `VirtualSpatialImage`.

## Value

SFE object with the new image added.

## Examples

``` r
library(EBImage)
#> 
#> Attaching package: ‘EBImage’
#> The following objects are masked from ‘package:SpatialFeatureExperiment’:
#> 
#>     affine, rotate, translate, transpose
library(SFEData)
library(RBioFormats)
#> BioFormats library version 7.3.0
fp <- tempfile()
fn <- XeniumOutput("v2", file_path = fp)
#> 
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> The downloaded files are in /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce82f5a1322/xenium2 
# Weirdly the first time I get the null pointer error
try(sfe <- readXenium(fn))
#> >>> Cell segmentations are found in `.parquet` file(s)
#> >>> Reading cell and nucleus segmentations
#> >>> Making MULTIPOLYGON nuclei geometries
#> >>> Making POLYGON cell geometries
#> Sanity checks on cell segmentation polygons:
#> >>> ..found 132 cells with (nested) polygon lists
#> >>> ..applying filtering
#> >>> Checking polygon validity
#> >>> Saving geometries to parquet files
#> >>> Reading cell metadata -> `cells.parquet`
#> >>> Reading h5 gene count matrix
#> >>> filtering cellSeg geometries to match 6272 cells with counts > 0
#> >>> filtering nucSeg geometries to match 6158 cells with counts > 0
sfe <- readXenium(fn)
#> >>> Preprocessed sf segmentations found
#> >>> Reading cell and nucleus segmentations
#> >>> Reading cell metadata -> `cells.parquet`
#> >>> Reading h5 gene count matrix
#> >>> filtering cellSeg geometries to match 6272 cells with counts > 0
#> >>> filtering nucSeg geometries to match 6158 cells with counts > 0
img <- getImg(sfe) |> toExtImage(resolution = 1L)
img <- img[,,1] > 500
Img(sfe, image_id = "mask") <- img
imageIDs(sfe)
#> [1] "morphology_focus" "mask"            
unlink(fn, recursive = TRUE)
```
