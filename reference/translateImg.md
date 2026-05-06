# Translate/shift image in space

This function shifts the spatial extent of the image in the x-y plane.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
translateImg(x, v, ...)

# S4 method for class 'BioFormatsImage'
translateImg(x, v, ...)

# S4 method for class 'ExtImage'
translateImg(x, v, ...)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- v:

  Numeric vector of length 2 to shift the image in the x-y plane.

- ...:

  Ignored. It's there so different methods can all be passed to the same
  `lapply` in the method for SFE objects. Some methods have extra
  arguments.

## Value

A `*Image` object of the same class that has been shifted in space.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
