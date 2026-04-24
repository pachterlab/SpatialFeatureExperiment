# Use the EBImage `Image` class in SFE objects

This is a thin wrapper around the
[`Image`](https://rdrr.io/pkg/EBImage/man/Image.html) class in the
`EBImage` package so it inherits from `VirtualSpatialImage` to be
compatible with `SpatialExperiment` from which SFE inherits. An `ext`
field is added to specify the spatial extent of the image in microns to
facilitate geometric operations on the SFE object (including the images)
and plotting with `Voyager`.

## Usage

``` r
# S4 method for class 'ExtImage'
show(object)

ExtImage(img, ext = NULL)
```

## Arguments

- object:

  An `ExtImage` object.

- img:

  An `Image` object or anything that inherits from `Image` such as
  `AnnotatedImage` in `RBioFormats`.

- ext:

  Numeric vector with names "xmin", "xmax", "ymin", "ymax" in microns
  indicating the spatial extent covered by the image. If `NULL`, then
  the extent will be inferred from the metadata, from physical pixel
  size and the number of pixels.

## Value

An `ExtImage` object.
