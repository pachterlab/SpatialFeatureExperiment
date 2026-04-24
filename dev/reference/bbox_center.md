# Find center of bounding box

Get x-y coordinates of the center of any bounding box

## Usage

``` r
bbox_center(bbox)
```

## Arguments

- bbox:

  A numeric vector of length 4 with names xmin, xmax, ymin, ymax, in any
  order.

## Value

A numeric vector of length 2.

## Examples

``` r
bbox <- c(xmin = 0, xmax = 100, ymin = 0, ymax = 80)
bbox_center(bbox)
#> [1] 50 40
```
