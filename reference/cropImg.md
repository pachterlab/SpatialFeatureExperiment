# Crop images

Crop images of class `*Image` in this package with a bounding box.

## Usage

``` r
# S4 method for class 'SpatRasterImage'
cropImg(x, bbox, filename = "")

# S4 method for class 'BioFormatsImage'
cropImg(x, bbox)

# S4 method for class 'ExtImage'
cropImg(x, bbox)
```

## Arguments

- x:

  An object of class `*Image` as implemented in this package.

- bbox:

  Numeric vector with names "xmin", "xmax", "ymin", "ymax", in any
  order, to specify the bounding box.

- filename:

  Output file name for transformed SpatRaster.

## Value

Image of the same class as input but cropped. For `BioFormatsImage`, the
image is not loaded into memory; only the extent is changed.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md),
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
