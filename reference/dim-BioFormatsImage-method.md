# Find dimension of BioFormatsImage

This is different from other classes. The metadata is read where the
dimensions in pixels can be found. The image itself is not read into
memory here.

## Usage

``` r
# S4 method for class 'BioFormatsImage'
dim(x)
```

## Arguments

- x:

  A
  [`BioFormatsImage`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/BioFormatsImage.md)
  object.

## Value

An integer vector of length 5 showing the number of rows and columns in
the full resolution image. The 5 dimensions are in the order of XYCZT:
x, y, channel, z, and time. This is not changed by transformations. Use
[`ext`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
to see the extent after transformation.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
