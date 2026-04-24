# Convert images to ExtImage

The `ExtImage` class is a thin wrapper around the `Image` class in
`ExtImage` so it inherits from `VirtualSpatialImage` as required by
`SpatialExperiment` and has extent as used in Voyager's plotting
functions. This function converts `SpatRasterImage` (thin wrapper around
the class in `terra`) and `BioFormatsImage` into `ExtImage` for image
operations as implemented in the `ExtImage` package.

## Usage

``` r
# S4 method for class 'BioFormatsImage'
toExtImage(x, resolution = 4L, channel = NULL)

# S4 method for class 'SpatRasterImage'
toExtImage(x, maxcell = 1e+07, channel = NULL)
```

## Arguments

- x:

  Either a `BioFormatsImage` or `SpatRasterImage` object.

- resolution:

  Integer, which resolution in the `BioFormatsImage` to read and
  convert. Defaults to 4, which is a lower resolution. Ignored if only 1
  resolution is present.

- channel:

  Integer vector to indicate channel(s) to read. If `NULL`, then all
  channels will be read.

- maxcell:

  Maximum number of pixels when `SpatRasterImage` is read into memory.

## Value

A `ExtImage` object. The image is loaded into memory.

## See also

toSpatRasterImage
