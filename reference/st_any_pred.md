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
# S4 method for class 'ANY'
st_any_pred(x, y, pred, yx = FALSE, sparse = FALSE, ...)

# S4 method for class 'duckspatial_df'
st_any_pred(x, y, pred, yx = NA, sparse = FALSE, ...)

# S4 method for class 'DuckDBTable'
st_any_pred(x, y, pred, yx = NA, sparse = TRUE, ...)

st_any_intersects(x, y, yx = FALSE, sparse = FALSE)

# S4 method for class 'ANY'
st_n_pred(x, y, pred, ...)

# S4 method for class 'duckspatial_df'
st_n_pred(x, y, pred, ...)

# S4 method for class 'DuckDBTable'
st_n_pred(x, y, pred, ...)

st_n_intersects(x, y)
```

## Arguments

- x:

  An object of class `sf`, `sfc`, or `sfg`.

- y:

  Another object of class `sf`, `sfc`, or `sfg`.

- pred:

  Name of a predicate such as "intersects", see
  [`ddbs_predicate`](https://cidree.github.io/duckspatial/reference/ddbs_predicate.html).
  It can correspond to functions in the `sf` package without the `st_`
  prefix. Deprecated behavior: A geometric binary predicate function,
  such as
  [`st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).

- yx:

  Whether to do `pred(y, x)` instead of `pred(x, y)`. For symmetric
  predicates, the results should be the same. When x has a large number
  of geometries and y has few, `pred(y, x)` is much faster than
  `pred(x, y)` for `st_intersects`, `st_disjoint`, and
  `st_is_within_distance`. Ignored for `DuckDBDataFrame`.

- sparse:

  If `TRUE`, returns numeric indices rather than logical vector.
  Defaults to `FALSE` for backward compatibility, though the default in
  `st_intersects` is `TRUE`.

- ...:

  Arguments passed to the predicate function, which would be `st_<pred>`
  for in memory geometries and
  [`ddbs_predicate`](https://cidree.github.io/duckspatial/reference/ddbs_predicate.html)
  when using DuckDB.

## Value

For `st_any_*`, a logical vector indicating whether each geometry in `x`
intersects (or other predicates such as is covered by) anything in `y`
or a numeric vector of indices of `TRUE` when `sparse = TRUE`.
Simplified from the `sgbp` results which indicate which item in `y` each
item in `x` intersects, which might not always be relevant. For
`st_n_*`, an integer vector indicating the number of geometries in y
returns TRUE for each geometry in x. For `DuckDB`, the output is a
`DuckDBDataColumn`.

## Examples

``` r
library(sf)
pts <- st_sfc(
    st_point(c(.5, .5)), st_point(c(1.5, 1.5)),
    st_point(c(2.5, 2.5))
)
pol <- st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
st_any_pred(pts, pol, pred = st_disjoint)
#> Warning: The `pred` argument of `st_any_pred()` should be a string as of
#> SpatialFeatureExperiment 1.16.0.
#> ℹ The deprecated feature was likely used in the SpatialFeatureExperiment
#>   package.
#>   Please report the issue at
#>   <https://github.com/pachterlab/SpatialFeatureExperiment/issues>.
#> [1] FALSE FALSE  TRUE
st_any_intersects(pts, pol)
#> [1]  TRUE  TRUE FALSE
st_n_pred(pts, pol, pred = st_disjoint)
#> [1] 0 0 1
st_n_intersects(pts, pol)
#> [1] 1 1 0
```
