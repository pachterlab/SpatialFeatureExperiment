# Scale image

This function scales the image about its center. After scaling, the
center of the image is not shifted.

## Usage

``` r
# S4 method for class 'AlignedSpatialImage'
scaleImg(x, factor, ...)

# S4 method for class 'BioFormatsImage'
scaleImg(x, factor, ...)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- factor:

  Numeric, scaling factor.

- ...:

  Ignored. It's there so different methods can all be passed to the same
  `lapply` in the method for SFE objects. Some methods have extra
  arguments.

## Value

A `*Image` object of the same class that has been scaled. Behind the
scene, it's only the extent that has been changed and the images are not
changed. The center of the image is unchanged.

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
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
