# Spatial graph methods

Spatial neighborhood graphs as `spdep`'s `listw` objects are stored in
the `int_metadata` of the SFE object. The `listw` class is used because
`spdep` has many useful methods that rely on the neighborhood graph as
`listw`. See the [`spdep` doumentation
website](https://r-spatial.github.io/spdep/reference/index.html) for
functions to edit the spatial neighborhood graph, or the `nb` object
within the `listw`.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
spatialGraphs(x, MARGIN = NULL, sample_id = "all", name = "all")

colGraphs(x, sample_id = "all", name = "all")

rowGraphs(x, sample_id = "all", name = "all")

annotGraphs(x, sample_id = "all", name = "all")

# S4 method for class 'SpatialFeatureExperiment'
spatialGraphs(x, MARGIN = NULL, sample_id = "all", name = "all") <- value

colGraphs(x, sample_id = "all", name = "all") <- value

rowGraphs(x, sample_id = "all", name = "all") <- value

annotGraphs(x, sample_id = "all", name = "all") <- value

# S4 method for class 'SpatialFeatureExperiment,numeric'
spatialGraphNames(x, MARGIN, sample_id = 1L)

# S4 method for class 'SpatialFeatureExperiment,numeric,ANY,character'
spatialGraphNames(x, MARGIN, sample_id = 1L) <- value

colGraphNames(x, sample_id = 1L)

rowGraphNames(x, sample_id = 1L)

annotGraphNames(x, sample_id = 1L)

colGraphNames(x, sample_id = 1L) <- value

rowGraphNames(x, sample_id = 1L) <- value

annotGraphNames(x, sample_id = 1L) <- value

# S4 method for class 'SpatialFeatureExperiment'
spatialGraph(x, type = 1L, MARGIN, sample_id = 1L)

colGraph(x, type = 1L, sample_id = 1L)

rowGraph(x, type = 1L, sample_id = 1L)

annotGraph(x, type = 1L, sample_id = 1L)

# S4 method for class 'SpatialFeatureExperiment'
spatialGraph(x, type = 1L, MARGIN, sample_id = NULL) <- value

colGraph(x, type = 1L, sample_id = 1L) <- value

rowGraph(x, type = 1L, sample_id = 1L) <- value

annotGraph(x, type = 1L, sample_id = 1L) <- value
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- MARGIN:

  As in [`apply`](https://rdrr.io/r/base/apply.html). 1 stands for rows
  and 2 stands for columns. In addition, 3 stands for spatial
  neighborhood graphs that correspond to `annotGeometries`.

- sample_id:

  Name of the sample the graph is associated with. This is useful when
  multiple pieces of tissues are in the same SFE object (say for a joint
  dimension reduction and clustering) and the spatial neighborhood is
  only meaningful within the same piece of tissue. See the `sample_id`
  argument in
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- name:

  Name of the graphs to add to each sample_id; used in the
  `spatialGraphs` replacement method as it must be character while
  `type` can be either an integer index or a name.

- value:

  A `listw` object (`*Graph`), or a named list of list of `listw`
  objects (`*Graphs`) where the names of the top level list are
  `sample_id`s when adding graphs for all samples in the margin of
  interest, or a list of `listw` objects when adding graphs for one
  sample in one margin.

- type:

  An integer specifying the index or string specifying the name of the
  \*Graph to query or replace. If missing, then the first item in the
  \*Graph will be returned or replaced.

## Value

Getters for multiple graphs return a named list. Getters for names
return a character vector of the names. Getters for single graphs return
a `listw` object. Setters return an SFE object.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
g1 <- findVisiumGraph(sfe)
g2 <- findSpatialNeighbors(sfe)

# Set all graphs of a margin by a named list
spatialGraphs(sfe, MARGIN = 2L, sample_id = "Vis5A") <-
    list(tri2nb = g2, visium = g1)
# Or equivalently
colGraphs(sfe, sample_id = "Vis5A") <- list(tri2nb = g2, visium = g1)

# Get all graphs of a margin, returning a named list
gs <- spatialGraphs(sfe, MARGIN = 2L)
# Or equivalently
gs <- colGraphs(sfe)

# Set graph of the same name and same margin for multiple samples
# Each sample has a separate graph
sfe2 <- McKellarMuscleData("small2")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe_combined <- cbind(sfe, sfe2)
colGraphs(sfe_combined, name = "visium", sample_id = "all") <-
    findVisiumGraph(sfe_combined, sample_id = "all")

# Get graph names
spatialGraphNames(sfe, MARGIN = 2L, sample_id = "Vis5A")
#> [1] "tri2nb" "visium"
# Or equivalently (sample_id optional as only one sample is present)
colGraphNames(sfe)
#> [1] "tri2nb" "visium"

# Set graph names
spatialGraphNames(sfe, MARGIN = 2L) <- c("foo", "bar")
colGraphNames(sfe) <- c("tri2nb", "visium")

# MARGIN = 1 means rowGraphs; MARGIN = 3 means annotation graphs (annotGraphs)
# for both getters and setters

# Set single graph by
# Spatial graph for myofibers
g_myofiber <- findSpatialNeighbors(sfe,
    type = "myofiber_simplified",
    MARGIN = 3L
)
spatialGraph(sfe, type = "myofiber", MARGIN = 3L) <- g_myofiber
# Or equivalently
annotGraph(sfe, "myofiber") <- g_myofiber

# Get a specific graph by name
g <- spatialGraph(sfe, "myofiber", MARGIN = 3L)
g2 <- spatialGraph(sfe, "visium", MARGIN = 2L)
# Or equivalently
g <- annotGraph(sfe, "myofiber")
g2 <- colGraph(sfe, "visium")
```
