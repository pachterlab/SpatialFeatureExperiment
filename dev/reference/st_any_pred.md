# Simple geometry predicates

Unlike functions in `sf` like `st_intersects`, this function simply
returns a logical vector indicating whether each geometry in `x`
intersects (or returns `TRUE` from other predicates) anything in `y`,
preferably when `y` only contains a small number of geometries or is one
single MULTI geometry. This is useful when cropping or subsetting an SFE
object with a geometry, such as tissue boundary or histological region
polygons or a bounding box.

## Usage

``` r
st_any_pred(x, y, pred, yx = FALSE, sparse = FALSE, ...)

st_any_intersects(x, y, yx = FALSE, sparse = FALSE)

st_n_pred(x, y, pred, ...)

st_n_intersects(x, y)
```

## Arguments

- x:

  An object of class `sf`, `sfc`, or `sfg`.

- y:

  Another object of class `sf`, `sfc`, or `sfg`.

- pred:

  A geometric binary predicate function, such as
  [`st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).
  It should return an object of class `sgbp`, for sparse predicates.

- yx:

  Whether to do `pred(y, x)` instead of `pred(x, y)`. For symmetric
  predicates, the results should be the same. When x has a large number
  of geometries and y has few, `pred(y, x)` is much faster than
  `pred(x, y)` for `st_intersects`, `st_disjoint`, and
  `st_is_within_distance`.

- sparse:

  If `TRUE`, returns numeric indices rather than logical vector.
  Defaults to `FALSE` for backward compatibility, though the default in
  `st_intersects` is `TRUE`.

- ...:

  Arguments passed to `pred`.

## Value

For `st_any_*`, a logical vector indicating whether each geometry in `x`
intersects (or other predicates such as is covered by) anything in `y`
or a numeric vector of indices of `TRUE` when `sparse = TRUE`.
Simplified from the `sgbp` results which indicate which item in `y` each
item in `x` intersects, which might not always be relevant. For
`st_n_*`, an integer vector indicating the number of geometries in y
returns TRUE for each geometry in x.

## Examples

``` r
library(sf)
#> Linking to GEOS 3.13.0, GDAL 3.8.5, PROJ 9.5.1; sf_use_s2() is TRUE
pts <- st_sfc(
    st_point(c(.5, .5)), st_point(c(1.5, 1.5)),
    st_point(c(2.5, 2.5))
)
pol <- st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
st_any_pred(pts, pol, pred = st_disjoint)
#> [1] FALSE FALSE  TRUE
st_any_intersects(pts, pol)
#> [1]  TRUE  TRUE FALSE
st_n_pred(pts, pol, pred = st_disjoint)
#> [1] 0 0 1
st_n_intersects(pts, pol)
#> [1] 1 1 0
```
