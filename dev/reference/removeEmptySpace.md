# Remove empty space

For each sample independently, all geometries and `spatialCoords` are
translated so the origin is at the minimum coordinates of the bounding
box of all geometries of the sample. This way coordinates of different
samples will be more comparable. This removes empty space in the images
if present.

## Usage

``` r
removeEmptySpace(sfe, sample_id = "all")
```

## Arguments

- sfe:

  An SFE object.

- sample_id:

  Sample to remove empty space.

## Value

An SFE object with empty space removed.

## Note

Unlike other functions in this package, this function operates on all
samples by default.

## Examples

``` r
library(SFEData)
library(SingleCellExperiment)
sfe <- McKellarMuscleData("full")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
# Only keep spots on tissue
sfe <- sfe[, colData(sfe)$in_tissue]
# Move the coordinates of the tissue
sfe <- removeEmptySpace(sfe)
```
