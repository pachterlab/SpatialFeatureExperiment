# Affine transformation of images

This function performs affine transformation on images, with any matrix
and translation vector.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
affineImg(x, M, v, maxcell = 1e+07, ...)

# S4 method for class 'BioFormatsImage'
affineImg(x, M, v, ...)

# S4 method for class 'ExtImage'
affineImg(x, M, v, ...)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- M:

  A 2x2 numeric matrix for the linear transformation in the xy plane.

- v:

  A numeric vector of length 2 for translation in the xy plane.

- maxcell:

  Max number of pixels to load `SpatRasterImage` into memory. The
  default 1e7 is chosen because this is the approximate number of pixels
  in the medium resolution image at `resolution = 4L` in Xenium OME-TIFF
  to make different methods of this function consistent.

- ...:

  Ignored. It's there so different methods can all be passed to the same
  `lapply` in the method for SFE objects. Some methods have extra
  arguments.

## Value

`SpatRasterImage` will be converted to `ExtImage`. Otherwise `*Image`
object of the same class. For `BioFormatsImage`, the transformation info
is stored and will be applied when the image is loaded into memory as
`ExtImage`.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
