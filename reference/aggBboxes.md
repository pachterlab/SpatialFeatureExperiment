# Aggregate bounding boxes

To find the bounding box of multiple bounding boxes.

## Usage

``` r
aggBboxes(bboxes)
```

## Arguments

- bboxes:

  Either a matrix with 4 rows whose columns are the different bounding
  boxes, with row names "xmin", "xmax", "ymin", and "ymax" in any order,
  or a list of bounding boxes which are named numeric vectors.

## Value

A named numeric vector for the total bounding box.

## Examples

``` r
bboxes <- list(c(xmin = 5, xmax = 10, ymin = 2, ymax = 20),
c(xmin = 8, xmax = 18, ymin = 0, ymax = 15))
bbox_all <- aggBboxes(bboxes)
```
