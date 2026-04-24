# Get the image from \*Image class

In SFE, S4 classes inheriting from `VirtualSpatialImage` have been
implemented to make these image classes compatible with
`SpatialExperiment`.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
imgRaster(x, maxcell = 1e+07, col = terra::map.pal("viridis", 100))

# S4 method for class 'BioFormatsImage'
imgRaster(x, resolution = 4L)

# S4 method for class 'ExtImage'
imgRaster(x)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- maxcell:

  positive integer. Maximum number of cells to use for the plot

- col:

  vector of colors. The default is `map.pal("viridis", 100)`

- resolution:

  Resolution to read in from OME-TIFF, defaults to 4, which is a medium
  resolution in Xenium.

## Value

Since version 1.9.8, `imgRaster` will return an array of hex colors, or
the raster object, as required by `SpatialExperiment`. This will break
older SFE code calling `imgRaster`.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ext.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/transposeImg.md)
