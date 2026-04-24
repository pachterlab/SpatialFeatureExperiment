# Source of images that are on disk

Get the file path of images that are on disk and not read into memory.
Only applies to `SpatRasterImage` and `BioFormatsImage`.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
imgSource(x)

# S4 method for class 'BioFormatsImage'
imgSource(x)

# S4 method for class 'ExtImage'
imgSource(x)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

## Value

String, file path to the original image on disk. For `SpatRasterImage`,
if the image is loaded into memory, then `NULL`.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/imgRaster.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/transposeImg.md)
