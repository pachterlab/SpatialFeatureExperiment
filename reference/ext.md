# Get and set extent of image objects

Unlike in `SpatialExperiment`, images in SFE have extents which are used
to align them to the geometries and in geometric operations on SFE
objects. These functions get or set the extent for S4 image classes
inheriting from `VirtualSpatialImage` implemented in the SFE package.

## Usage

``` r
# S4 method for class 'BioFormatsImage'
ext(x)

# S4 method for class 'ExtImage'
ext(x)

# S4 method for class 'SpatRasterImage'
ext(x)

# S4 method for class 'BioFormatsImage,numeric'
ext(x) <- value

# S4 method for class 'ExtImage,numeric'
ext(x) <- value

# S4 method for class 'SpatRasterImage,numeric'
ext(x) <- value
```

## Arguments

- x:

  A `*Image` object.

- value:

  A numeric vector with names "xmin", "xmax", "ymin", "ymax" specifying
  the extent to use.

## Value

Getters return a numeric vector specifying the extent. Setters return a
`*Image` object of the same class as the input.

## Note

For `SpatRasterImage`, the image may be may not be loaded into memory.
You can check if the image is loaded into memory with
`terra::inMemory(x)`, and check the original file path with
[`imgSource`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md).
If the image is not loaded into memory, then the original file must be
present at the path indicated by
[`imgSource`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md)
in order for any code using the image to work, which includes this
function `ext`.

For `BioFormatsImage`, internally only the pre-transform extent is
stored. The `ext` getter will apply the transformation on the fly. The
setter sets the pre-transformation extent.

## See also

Other image methods:
[`SFE-image`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md),
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-ExtImage-method.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
