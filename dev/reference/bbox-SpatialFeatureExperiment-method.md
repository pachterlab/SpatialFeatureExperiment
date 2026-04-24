# Find bounding box of SFE objects

Find bounding box of the union of all `colGeometries` and
`annotGeometries` of each sample in the SFE object. This can be used to
remove empty space so the tissue and geometries have one corner at the
origin so all samples will be on comparable coordinates.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
bbox(sfe, sample_id = "all", include_images = FALSE, include_row = TRUE)
```

## Arguments

- sfe:

  A `SpatialFeatureExperiment` object.

- sample_id:

  Sample(s) whose bounding box(es) to find. The bounding box would be
  for the union of all `colGeometries` and `annotGeometries` associated
  with each sample.

- include_images:

  Logical, whether the bounding boxes should include image extents.
  Defaults to `FALSE` because often the image has a lot of empty space
  surrounding the tissue.

- include_row:

  Logical, whether the bounding boxes should include `rowGeometries`,
  defaults to `TRUE`.

## Value

For one sample, then a named vector with names `xmin`, `ymin`, `xmax`,
and `ymax` specifying the bounding box. For multiple samples, then a
matrix whose columns are samples and whose rows delineate the bounding
box.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
bbox(sfe, sample_id = "Vis5A")
#>  xmin  ymin  xmax  ymax 
#>  5000 13000  7000 15000 
```
