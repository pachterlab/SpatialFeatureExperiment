# Print method for SpatialFeatureExperiment

Printing summaries of `colGeometries`, `rowGeometries`, and
`annotGeometries` in addition to what's shown for `SpatialExperiment`.
Geometry names and types are printed.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
show(object)
```

## Arguments

- object:

  A `SpatialFeatureExperiment` object.

## Value

None (invisible `NULL`).

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe # The show method is implicitly called
#> class: SpatialFeatureExperiment 
#> dim: 15123 77 
#> metadata(0):
#> assays(1): counts
#> rownames(15123): ENSMUSG00000025902 ENSMUSG00000096126 ...
#>   ENSMUSG00000064368 ENSMUSG00000064370
#> rowData names(6): Ensembl symbol ... vars cv2
#> colnames(77): AAATTACCTATCGATG AACATATCAACTGGTG ... TTCTTTGGTCGCGACG
#>   TTGATGTGTAGTCCCG
#> colData names(12): barcode col ... prop_mito in_tissue
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : imageX imageY
#> imgData names(1): sample_id
#> 
#> unit: full_res_image_pixels
#> Geometries:
#> colGeometries: spotPoly (POLYGON) 
#> annotGeometries: tissueBoundary (POLYGON), myofiber_full (GEOMETRY), myofiber_simplified (GEOMETRY), nuclei (POLYGON), nuclei_centroid (POINT) 
#> 
#> Graphs:
#> Vis5A: 
```
