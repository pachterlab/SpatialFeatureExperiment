# SpatRaster representation of images in SFE objects

`SpatialFeatureExperiment` and the `Voyager` package work with images
differently from `SpatialExperiment`. In SFE and `Voyager`'s, plotting
functions for SFE objects, the images can be read with
[`rast`](https://rspatial.github.io/terra/reference/rast.html) and
represented as `SpatRaster`, so the image is not entirely loaded into
memory unless necessary. Plotting will not load a large image into
memory; rather the image will be downsampled and the downsampled version
is plotted. A `SpatRasterImage` object (as of Bioc 3.19 or SFE version
1.6 and above) is a `SpatRaster` object but also inheriting from
`VirtualSpatialImage` as required by `SpatialExperiment`.

## Usage

``` r
SpatRasterImage(img)

# S4 method for class 'SpatRasterImage'
show(object)
```

## Arguments

- img:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or `PackedSpatRaster` object.

- object:

  A `SpatRasterImage` object.

## Value

A `SpatRasterImage` object.

## Examples

``` r
# Example code
```
