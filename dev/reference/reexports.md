# Functions re-exported from other packages

These are some commonly used getters and setters of classes that SFE
inherits so you don't have to separately attach those packages to use
these functions.

## Usage

``` r
colData(x, ...)

rowData(x, use.names = TRUE, ...)

colData(x, ...) <- value

spatialCoords(x, ...)

spatialCoords(x) <- value

spatialCoordsNames(x)

getImg(x, ...)

imgData(x)

rmvImg(x, ...)

counts(object, ...)

logcounts(object, ...)

reducedDim(x, type, ...)
```

## Arguments

- x:

  A SummarizedExperiment object or derivative.

- ...:

  For `assay`, arguments in `...` are forwarded to `assays`.

  For `rbind`, `cbind`, `...` contains SummarizedExperiment objects (or
  derivatives) to be combined.

  For other accessors, ignored.

- use.names:

  For `rowData`: Like
  [`mcols`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html)`(x)`,
  by default `rowData(x)` propagates the rownames of `x` to the returned
  [DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  object (note that for a SummarizedExperiment object or derivative, the
  rownames are also the names i.e. `rownames(x)` is always the same as
  `names(x)`). Setting `use.names=FALSE` suppresses this propagation
  i.e. it returns a
  [DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  object with no rownames. Use this when `rowData(x)` fails, which can
  happen when the rownames contain NAs (because the rownames of a
  SummarizedExperiment object or derivative can contain NAs, but the
  rownames of a
  [DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  object cannot).

  For `combineRows` and `combineCols`: See Combining section below.

- value:

  An object of a class specified in the S4 method signature or as
  outlined in ‘Details’.

- object:

  A `SingleCellExperiment` object, which includes SFE.

- type:

  Name or numeric index to indicate which `reducedDim` to get, such as
  "PCA". By default the first item in `reducedDims`.
