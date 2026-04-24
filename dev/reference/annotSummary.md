# Summarize attributes of an annotGeometry for each cell/spot

In SFE objects, the annotation geometries don't have to correspond to
the dimensions of the gene count matrix, so there generally is no one to
one mapping between annotation geometries and cells/spots. However, it
may be interesting to relate attributes of annotation geometries to
cell/spots so the attributes can be related to gene expression. This
function summarizes attributes of an `annotGeometry` for each cell/spot
by a geometric predicate with a `colGeometry`.

## Usage

``` r
annotSummary(
  sfe,
  colGeometryName = 1L,
  annotGeometryName = 1L,
  annotColNames = 1L,
  sample_id = "all",
  pred = st_intersects,
  summary_fun = mean
)
```

## Arguments

- sfe:

  An SFE object.

- colGeometryName:

  Name of column geometry for the predicate.

- annotGeometryName:

  Name of annotation geometry for the predicate.

- annotColNames:

  Character, column names of the `annotGeometry` of interest, to
  indicate the columns to summarize. Columns that are absent from the
  `annotGeometry` are removed. The column cannot be "geometry" or
  "barcode".

- sample_id:

  Which sample(s) to operate on. Can be "all" to indicate all samples.

- pred:

  Predicate function to use, defaults to
  [`st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html).

- summary_fun:

  Function for the summary, defaults to `mean`.

## Value

A data frame whose row names are the relevant column names of `sfe`, and
each column of which is the summary of each column specified in
`annotColName`.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
s <- annotSummary(sfe, "spotPoly", "myofiber_simplified",
    annotColNames = c("area", "convexity")
)
```
