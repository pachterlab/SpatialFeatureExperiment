# Find dimensions of ExtImage

This method exists to make the output of
[`dim()`](https://rdrr.io/r/base/dim.html) for `ExtImage` consistent
with that of `Image` which `ExtImage` inherits from, overriding the
`VirtualSpatialImage` method.

## Usage

``` r
# S4 method for class 'ExtImage'
dim(x)
```

## Arguments

- x:

  A
  [`ExtImage`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ExtImage.md)
  object.

## Value

An integer vector. As in `EBImage`, the first element indicates number
of pixels in the x direction, or number of columns in the image, and the
second element indicates the number of pixels in the y direction. This
is unlike array indexing.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-BioFormatsImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
