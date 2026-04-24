# Mirror/flip images

Flip images along the middle horizontal or vertical axis.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
mirrorImg(
  x,
  direction = c("vertical", "horizontal"),
  filename = "",
  maxcell = NULL,
  ...
)

# S4 method for class 'BioFormatsImage'
mirrorImg(x, direction = c("vertical", "horizontal"), ...)

# S4 method for class 'ExtImage'
mirrorImg(x, direction = c("vertical", "horizontal"), ...)
```

## Arguments

- x:

  SpatRaster or SpatVector

- direction:

  character. Should (partially) match "vertical" to flip by rows, or
  "horizontal" to flip by columns

- filename:

  character. Output filename

- maxcell:

  Max number of pixels to load `SpatRasterImage` into memory. The
  default 1e7 is chosen because this is the approximate number of pixels
  in the medium resolution image at `resolution = 4L` in Xenium OME-TIFF
  to make different methods of this function consistent.

- ...:

  additional arguments for writing files as in
  [`writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.html)

## Value

`*Image` object of the same class.

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
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/transposeImg.md)
