# Get and set results from local spatial statistics

Local spatial statics like local Moran's I, local Geary's C, Getis-Ord
Gi\*, and geographically weighted summary statistics return values at
each spatial location. Just like dimension reductions, these results are
clearly associated with the broader SFE object, so they should have a
place within the object. However, a separate field is needed because
these analyses are conceptually distinct from dimension reduction. Also,
each feature (e.g. gene) can have its own results with values at each
location. The `localResults` field in the SFE object stores these
results that has a value for each spatial location.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
localResults(
  x,
  sample_id = "all",
  name = "all",
  features = NULL,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  withDimnames = TRUE,
  swap_rownames = NULL,
  ...
)

# S4 method for class 'SpatialFeatureExperiment'
localResults(
  x,
  sample_id = "all",
  name = "all",
  features = NULL,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  withDimnames = TRUE,
  swap_rownames = NULL,
  ...
) <- value

# S4 method for class 'SpatialFeatureExperiment'
localResultNames(x)

# S4 method for class 'SpatialFeatureExperiment,character'
localResultNames(x) <- value

# S4 method for class 'SpatialFeatureExperiment'
localResultFeatures(
  x,
  type = 1L,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  swap_rownames = NULL
)

# S4 method for class 'SpatialFeatureExperiment'
localResultAttrs(
  x,
  type = 1L,
  feature,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  swap_rownames = NULL
)

# S4 method for class 'SpatialFeatureExperiment'
localResult(
  x,
  type = 1L,
  feature,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  sample_id = 1L,
  withDimnames = TRUE,
  simplify = TRUE,
  swap_rownames = NULL
)

# S4 method for class 'SpatialFeatureExperiment'
localResult(
  x,
  type = 1L,
  feature,
  colGeometryName = NULL,
  annotGeometryName = NULL,
  sample_id = 1L,
  withDimnames = TRUE
) <- value
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- sample_id:

  Sample ID to get or set geometries.

- name:

  Name of the spatial method used, such as "localmoran".

- features:

  Features whose local results to get or set, for `localResults` getter
  and setter for multiple features at a time.

- colGeometryName:

  Which `colGeometry` to get or set local results.

- annotGeometryName:

  Which `annotGeometry` to get or set local results.

- withDimnames:

  Logical. If `TRUE`, then the dimnames (colnames or rownames) of the
  gene count matrix should correspond to row names of the `sf` data
  frames of interest.

- swap_rownames:

  Name of a column in `rowData` to identify features instead of the row
  names of the SFE object. For example, if the row names of the SFE
  object are Ensembl IDs and gene symbols are in the "symbol" column in
  `rowData`, then putting "symbol" for this argument will use the gene
  symbols to identify which gene's local results to get or set.

- ...:

  Ignored

- value:

  Values to set, should be either a matrix or a data frame.

- type:

  Name or index of the spatial method used, such as "localmoran".

- feature:

  Feature whose local results to get or set, for `localResult` getter
  and setter for one feature at a time.

- simplify:

  Basically whether to return the content of the list rather than a list
  when the list only has one element, such as results for one type and
  one feature.

## Value

`localResults` returns a named list each element of which is a set of
local results of interest. `localResult` returns a matrix or a data
frame, whichever the original is when it's set. `localResultNames`
returns a character vector. Setters return an SFE object with the
desired field set. For genes and `colData` columns, the local results
are stored in the `localResults` field in `int_colData`, whereas for
`colGeometries` and `annotGeometries`, the local results are stored as
columns in the same `sf` data frames. `localResultFeatures` returns a
character vector of names of features for which local results are
available. `localResultAttrs` returns a character vector of the column
names of the local results of one type for one feature. It returns
`NULL` if the results are a vector.

## Examples

``` r
# Toy example
sfe <- readRDS(system.file("extdata/sfe_toy.rds",
    package = "SpatialFeatureExperiment"
))
# localResults functions are written for organizing results from local
# spatial statistics (see the Voyager package). But for the examples here,
# random toy matrices are used. The real results are often matrices, with a
# matrix for each feature.
library(S4Vectors)
set.seed(29)
toy_res1 <- matrix(rnorm(10),
    nrow = 5, ncol = 2,
    dimnames = list(colnames(sfe), c("meow", "purr"))
)
toy_res1b <- matrix(rgamma(10, shape = 2),
    nrow = 5, ncol = 2,
    dimnames = list(colnames(sfe), c("meow", "purr"))
)
toy_df1 <- DataFrame(gene1 = I(toy_res1), gene2 = I(toy_res1b))

toy_res2 <- matrix(rpois(10, lambda = 2),
    nrow = 5, ncol = 2,
    dimnames = list(colnames(sfe), c("sassy", "tortitude"))
)
toy_df2 <- DataFrame(gene1 = I(toy_res2))
# Set all local results
localResults(sfe) <- list(localmoran = toy_df1, Gistar = toy_df2)
# Get all local results
lrs <- localResults(sfe)

# Set results of the same type for multiple genes
localResults(sfe, name = "localmoran") <- toy_df1
# Can also use a list
localResults(sfe, name = "localmoran") <- as.list(toy_df1)
# Get results of the same type for multiple genes
lrs <- localResults(sfe, name = "localmoran", features = c("gene1", "gene2"))

# Set results for one type and one gene
localResult(sfe, "localmoran", feature = "gene1") <- toy_res1
# Get results for one type and one gene
lr <- localResult(sfe, "localmoran", feature = "gene1")

# Set results for a feature in colGeometries
cg_toy <- readRDS(system.file("extdata/cg_toy.rds",
    package = "SpatialFeatureExperiment"
))
colGeometry(sfe, "cg") <- cg_toy
localResult(sfe, "localmoran",
    feature = "gene1",
    colGeometryName = "cg"
) <- toy_res1
# Get results for a feature in colGeometries
lr <- localResult(sfe, "localmoran", "gene1", colGeometryName = "cg")
```
