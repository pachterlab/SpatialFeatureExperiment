# Convert multiple listw graphs into a single sparse adjacency matrix

Each sample in the SFE object has a separate spatial neighborhood graph.
Spatial analyses performed jointly on multiple samples require a
combined spatial neighborhood graph from the different samples, where
the different samples would be disconnected components of the graph.
This combined adjacency matrix can be used in MULTISPATI PCA.

## Usage

``` r
multi_listw2sparse(listws)
```

## Arguments

- listws:

  A list of `listw` objects.

## Value

A sparse `dgCMatrix` of the combined spatial neighborhood graph, with
the original spatial neighborhood graphs of the samples on the diagonal.
When the input is an SFE object, the rows and columns will match the
column names of the SFE object.
