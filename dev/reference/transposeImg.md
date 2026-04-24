# Transpose images

Swap rows and columns of images. In effect, this will flip the image
around the diagonal running from top left to bottom right.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
transposeImg(x, filename = "", maxcell = NULL, ...)

# S4 method for class 'BioFormatsImage'
transposeImg(x, ...)

# S4 method for class 'ExtImage'
transposeImg(x, ...)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- filename:

  Output file name for transformed SpatRaster.

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

For `SpatRasterImage` and `ExtImage`, object of the same class. For
`BioFormatsImage`, the image of the specified resolution is read into
memory and then the `ExtImage` method is called, returning `ExtImage`.
For the extent: xmin and xmax are switched with ymin and ymax.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/translateImg.md)
