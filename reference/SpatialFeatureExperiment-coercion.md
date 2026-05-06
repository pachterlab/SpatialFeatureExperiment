# SpatialFeatureExperiment coercion methods

The `SpatialFeatureExperiment` class inherits from `SpatialExperiment`,
which in turn inherits from `SingleCellExperiment`. A
`SpatialExperiment` object with geometries in `colGeometries` in the
`int_colData`, `rowGeometries` in the `int_elementMetadata`, or
`annotGeometries` in the `int_metadata` can be directly converted to
`SpatialFeatureExperiment` with `as(spe, "SpatialFeatureExperiment")`. A
`SpatialExperiment` object without the geometries can also be converted;
the coordinates in the `spatialCoords` field will be used to make POINT
geometries named "centroids" to add to `colGeometries`. The geometries
can also be supplied separately when using `toSpatialFeatureExperiment`.
Images are converted to `SpatRaster`.

## Usage

``` r
# S4 method for class 'SpatialExperiment'
toSpatialFeatureExperiment(
  x,
  colGeometries = NULL,
  rowGeometries = NULL,
  annotGeometries = NULL,
  spatialCoordsNames = c("x", "y"),
  annotGeometryType = "POLYGON",
  spatialGraphs = NULL,
  spotDiameter = NA,
  unit = NULL
)

# S4 method for class 'SingleCellExperiment'
toSpatialFeatureExperiment(
  x,
  sample_id = "sample01",
  spatialCoordsNames = c("x", "y"),
  spatialCoords = NULL,
  colGeometries = NULL,
  rowGeometries = NULL,
  annotGeometries = NULL,
  annotGeometryType = "POLYGON",
  spatialGraphs = NULL,
  spotDiameter = NA,
  scaleFactors = 1,
  imageSources = NULL,
  image_id = NULL,
  loadImage = TRUE,
  imgData = NULL,
  unit = NULL
)

# S4 method for class 'Seurat'
toSpatialFeatureExperiment(
  x,
  add_molecules = TRUE,
  flip = c("geometry", "image", "none"),
  image_scalefactors = c("lowres", "hires"),
  unit = NULL,
  BPPARAM = SerialParam()
)
```

## Arguments

- x:

  A `SpatialExperiment` or `Seurat` object to be coerced to a
  `SpatialFeatureExperiment` object.

- colGeometries:

  Geometry of the entities that correspond to the columns of the gene
  count matrix, such as cells and Visium spots. It must be a named list
  of one of the following:

  An `sf` data frame

  :   The geometry column specifies the geometry of the entities.

  An ordinary data frame specifying centroids

  :   Column names for the coordinates are specified in the
      `spatialCoordsNames` argument. For Visium and ST, in addition to
      the centroid coordinate data frame, the spot diameter in the same
      unit as the coordinates can be specified in the `spotDiamter`
      argument.

  An ordinary data frame specifying polygons

  :   Also use `spatialCoordsNames`. There should an additional column
      "ID" to specify which vertices belong to which polygon. The
      coordinates should not be in list columns. Rather, the data frame
      should look like it is passed to
      [`ggplot2::geom_polygon`](https://ggplot2.tidyverse.org/reference/geom_polygon.html).
      If there are holes, then there must also be a column "subID" that
      differentiates between the outer polygon and the holes.

  In all cases, the data frame should specify the same number of
  geometries as the number of columns in the gene count matrix. If the
  column "barcode" is present, then it will be matched to column names
  of the gene count matrix. Otherwise, the geometries are assumed to be
  in the same order as columns in the gene count matrix. If the
  geometries are specified in an ordinary data frame, then it will be
  converted into `sf` internally. Named list of data frames because each
  entity can have multiple geometries, such as whole cell and nuclei
  segmentations. The geometries are assumed to be POINTs for centroids
  and POLYGONs for segmentations. If polygons are specified in an
  ordinary data frame, then anything with fewer than 3 vertices will be
  removed. For anything other than POINTs, attributes of the geometry
  will be ignored.

- rowGeometries:

  Geometry associated with genes or features, which correspond to rows
  of the gene count matrix.

- annotGeometries:

  Geometry of entities that do not correspond to columns or rows of the
  gene count matrix, such as tissue boundary and pathologist annotations
  of histological regions, and nuclei segmentation in a Visium dataset.
  Also a named list as in `colGeometries`. The ordinary data frame may
  specify POINTs, POLYGONs, or LINESTRINGs, or their MULTI versions.
  Each data frame can only specify one type of geometry. For MULTI
  versions, there must be a column "group" to identify each MULTI
  geometry.

- spatialCoordsNames:

  A `character` vector of column names if `*Geometries` arguments have
  ordinary data frames, to identify the columns in the ordinary data
  frames that specify the spatial coordinates. If `colGeometries` is not
  specified, then this argument will behave as in
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html),
  but `colGeometries` will be given precedence if provided.

- annotGeometryType:

  Character vector specifying geometry type of each element of the list
  if `annotGeometry` is specified. Each element of the vector must be
  one of POINT, LINESTRING, POLYGON, MULTIPOINT, MULTILINESTRING, and
  MULTIPOLYGON. Must be either length 1 (same for all elements of the
  list) or the same length as the list. Ignored if the corresponding
  element is an `sf` object.

- spatialGraphs:

  A named list of `listw` objects (see `spdep`) for spatial neighborhood
  graphs.

- spotDiameter:

  Spot diameter for technologies with arrays of spots of fixed diameter
  per slide, such as Visium, ST, DBiT-seq, and slide-seq. The diameter
  must be in the same unit as the coordinates in the \*Geometry
  arguments. Ignored for geometries that are not POINT or MULTIPOINT.

- unit:

  \# Default unit is `"micron"`. However for Visium one can choose
  between `"micron"` or `"full_res_image_pixel"`.

- sample_id:

  A `character` sample identifier, which matches the `sample_id` in
  [`imgData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html).
  The `sample_id` will also be stored in a new column in
  [`colData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-colData.html),
  if not already present. Default = `sample01`.

- spatialCoords:

  A numeric matrix containing columns of spatial coordinates, as in
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).
  The coordinates are centroids of the entities represented by the
  columns of the gene count matrix. If `colGeometries` is also
  specified, then it will be given priority and a warning is issued.
  Otherwise, the `sf` representation of the centroids will be stored in
  the `colGeometry` called `centroids`.

- scaleFactors:

  Optional scale factors associated with the image(s). This can be
  provided as a numeric value, numeric vector, list, or file path to a
  JSON file for the 10x Genomics Visium platform. For 10x Genomics
  Visium, the correct scale factor will automatically be selected
  depending on the resolution of the image from `imageSources`. Default
  = `1`.

- imageSources:

  Optional file path(s) or URL(s) for one or more image sources.

- image_id:

  Optional character vector (same length as `imageSources`) containing
  unique image identifiers.

- loadImage:

  Logical indicating whether to load image into memory. Default =
  `FALSE`.

- imgData:

  Optional `DataFrame` containing the image data. Alternatively, this
  can be built from the arguments `imageSources` and `image_id` (see
  Details).

- add_molecules:

  Logical, whether to add transcripts coordinates to an object.

- flip:

  To flip the image, geometry coordinates, or none. Because the image
  has the origin at the top left while the geometry has origin at the
  bottom left, one of them needs to be flipped for them to match. If one
  of them is already flipped, then use "none". The image will not be
  flipped if it's GeoTIFF.

- image_scalefactors:

  \# A `character`, choose between "lowres" or "hires". Only for 10X
  Visium, image scaling factors are from \`scalefactors_json.json\`
  file.

- BPPARAM:

  Deprecated when coercing from `SpatialExperiment`, but is used when
  coercing from `Seurat` object.

## Value

A `SpatialFeatureExperiment` object

## Examples

``` r
library(VisiumIO)
#> Loading required package: TENxIO
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> 
#> Attaching package: ‘IRanges’
#> The following objects are masked from ‘package:EBImage’:
#> 
#>     resize, tile
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
#> The following object is masked from ‘package:EBImage’:
#> 
#>     channel
# From examples of TENxVisium()
sample_dir <- system.file(
file.path("extdata", "10xVisium", "section1"),
package = "VisiumIO"
)
## using spacerangerOut folder
tv <- TENxVisium(
    spacerangerOut = file.path(sample_dir, "outs"), processing = "raw", images = "lowres"
)
spe <- import(tv)
# There can't be duplicate barcodes
colnames(spe) <- make.unique(colnames(spe), sep = "-")
rownames(spatialCoords(spe)) <- colnames(spe)
sfe <- toSpatialFeatureExperiment(spe)
# For coercing Seurat to SFE see this -> ./vignettes/seurat_sfe_coerce.Rmd
```
