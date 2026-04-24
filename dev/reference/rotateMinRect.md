# Rotate the tissue to reduce empty space in plots

This function finds the minimum bounding rotated rectangle and then
rotates the sample so the minimum bounding rotated rectangle is
horizontal or vertical. This will reduce empty space in plots when the
sample has a long shape.

## Usage

``` r
# S4 method for class 'matrix'
rotateMinRect(x, orientation = c("horizontal", "vertical"), ...)

# S4 method for class 'sfc'
rotateMinRect(x, orientation = c("horizontal", "vertical"), ...)

# S4 method for class 'sf'
rotateMinRect(x, orientation = c("horizontal", "vertical"), ...)

# S4 method for class 'SpatialFeatureExperiment'
rotateMinRect(x, orientation = c("horizontal", "vertical"), maxcell = 1e+07)
```

## Arguments

- x:

  Object of interest to be rotated, can be a matrix with 2 columns,
  `sf`, `sfc`, or `SpatialFeatureExperiment`.

- orientation:

  Whether you want the rotated sample to be horizontal or vertical.

- ...:

  Ignored

- maxcell:

  Rotating `SpatRasterImage` will convert it into `ExtImage`, loading
  the image into memory. This argument specifies the maximum number of
  pixels in the image loaded into memory. The image will be down sampled
  to approximately this number of pixels.

## Value

Object of the same class as `x`, rotated.
