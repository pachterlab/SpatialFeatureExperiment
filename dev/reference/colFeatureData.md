# Get global spatial analysis results and metadata of colData, rowData, and geometries

Results of spatial analyses on columns in `colData`, `rowData`, and
geometries are stored in their metadata. The `colFeaturedata` function
allows the users to more directly access these results.

## Usage

``` r
colFeatureData(sfe)

rowFeatureData(sfe)

geometryFeatureData(sfe, type, MARGIN = 2L)

reducedDimFeatureData(sfe, dimred)
```

## Arguments

- sfe:

  An SFE object.

- type:

  Which geometry, can be name (character) or index (integer)

- MARGIN:

  Integer, 1 means rowGeometry, 2 means colGeometry, and 3 means
  annotGeometry. Defaults to 2, colGeometry.

- dimred:

  Name of a dimension reduction, can be seen in
  [`reducedDimNames`](https://rdrr.io/pkg/SingleCellExperiment/man/reducedDims.html).

## Value

A `DataFrame`.

## See also

getParams

## Examples

``` r
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SFEData)
library(Voyager)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
colGraph(sfe, "visium") <- findVisiumGraph(sfe)
# Moran's I for colData
sfe <- colDataMoransI(sfe, "nCounts")
colFeatureData(sfe)
#> DataFrame with 12 rows and 2 columns
#>           moran_Vis5A   K_Vis5A
#>             <numeric> <numeric>
#> barcode            NA        NA
#> col                NA        NA
#> row                NA        NA
#> x                  NA        NA
#> y                  NA        NA
#> ...               ...       ...
#> sample_id          NA        NA
#> nCounts      0.675416   1.67027
#> nGenes             NA        NA
#> prop_mito          NA        NA
#> in_tissue          NA        NA
```
