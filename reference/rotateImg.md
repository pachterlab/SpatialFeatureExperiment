# Rotate image

As in `SpatialExperiment`, rotation here must be a multiple of 90
degrees.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
rotateImg(x, degrees, maxcell = 1e+07, ...)

# S4 method for class 'BioFormatsImage'
rotateImg(x, degrees, ...)

# S4 method for class 'ExtImage'
rotateImg(x, degrees, ...)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- degrees:

  How many degrees to rotate. Positive number means clockwise and
  negative number means counterclockwise.

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

`SpatRasterImage` will be loaded into memory and converted to
`ExtImage`. Otherwise `*Image` object of the same class.

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
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
