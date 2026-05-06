# Row geometry getters and setters

`rowGeometries` are geometries that corresponding to rows of the gene
count matrix, such as smFISH transcript spots. The `txSpots()` function
is a convenience wrapper for transcript spots, although this entirely
depends on the `rowGeometry` being named `txSpots`.

## Usage

``` r
rowGeometry(x, type = 1L, sample_id = 1L, withDimnames = TRUE)

rowGeometry(
  x,
  type = 1L,
  sample_id = 1L,
  withDimnames = TRUE,
  partial = FALSE,
  translate = TRUE
) <- value

rowGeometries(x, sample_id = "all", withDimnames = TRUE)

rowGeometries(
  x,
  sample_id = "all",
  withDimnames = TRUE,
  partial = FALSE,
  translate = TRUE
) <- value

rowGeometryNames(x)

rowGeometryNames(x) <- value

txSpots(x, sample_id = 1L, withDimnames = TRUE)

txSpots(
  x,
  sample_id = 1L,
  withDimnames = TRUE,
  partial = FALSE,
  translate = TRUE
) <- value
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- type:

  An integer specifying the index or string specifying the name of the
  \*Geometry to query or replace. If missing, then the first item in the
  \*Geometries will be returned or replaced.

- sample_id:

  Sample ID to get or set geometries.

- withDimnames:

  Logical. If `TRUE`, then the dimnames (colnames or rownames) of the
  gene count matrix should correspond to row names of the `sf` data
  frames of interest.

- partial:

  In setters, if a `rowGeometry` of the same name exists, whether to
  only replace the rows present in `value`.

- translate:

  Logical. Only used if
  [`removeEmptySpace`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/removeEmptySpace.md)
  has been run of the SFE object. If that's the case, this argument
  indicates whether the new value to be assigned to the geometry is in
  the coordinates prior to removal of empty space so it should be
  translated to match the new coordinates after removing empty space.
  Default to `TRUE`.

- value:

  Value to set. For `dimGeometry`, must be a `sf` data frame with the
  same number of rows as size in the dimension of interest, or an
  ordinary data frame that can be converted to such a `sf` data frame
  (see
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md)).
  For `dimGeometries`, must be a list of such `sf` or ordinary data
  frames.

## Details

When there are multiple samples in the SFE object, `rowGeometries` for
each sample has the `sample_id` appended to the name of the geometry.
For example, if the name is `txSpots` and the sample ID is `sample01`,
then the actual name of the `rowGeometry` is `txSpots_sample01`. In the
getter, one can still specify
`rowGeometry(sfe, "txSpots", sample_id = "sample01")`.

Appending the `sample_id` is unnecessary when there is only one sample,
but `sample_id` will be appended when to SFE objects are combined with
`cbind`. It is necessary to distinguish bewteen different samples
because they can have overlapping coordinate values.

## See also

\[dimGeometries()\], \[colGeometries()\]

## Examples

``` r
library(SFEData)
library(RBioFormats)
fp <- tempfile()
dir_use <- XeniumOutput("v2", file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> The downloaded files are in /tmp/RtmphTOllm/fileea9d3dabca41/xenium2 
# RBioFormats issue
try(sfe <- readXenium(dir_use, add_molecules = TRUE))
#> >>> Must use gene symbols as row names when adding transcript spots.
#> >>> Cell segmentations are found in `.parquet` file(s)
#> >>> Reading cell and nucleus segmentations
#> >>> Making MULTIPOLYGON nuclei geometries
#> >>> Making POLYGON cell geometries
#> Sanity checks on cell segmentation polygons:
#> >>> ..found 132 cells with (nested) polygon lists
#> >>> ..applying filtering
#> >>> Checking polygon validity
#> >>> Saving geometries to parquet files
#> >>> Reading cell metadata -> `cells.parquet`
#> >>> Reading h5 gene count matrix
#> >>> filtering cellSeg geometries to match 6272 cells with counts > 0
#> >>> filtering nucSeg geometries to match 6158 cells with counts > 0
#> >>> Reading transcript coordinates
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk
#> >>> Total of 116 features/genes with no transcript detected or `min_phred` < 20 are removed from SFE object
#> >>> To keep all features -> set `min_phred = NULL`
sfe <- readXenium(dir_use, add_molecules = TRUE)
#> >>> Must use gene symbols as row names when adding transcript spots.
#> >>> Preprocessed sf segmentations found
#> >>> Reading cell and nucleus segmentations
#> >>> Reading cell metadata -> `cells.parquet`
#> >>> Reading h5 gene count matrix
#> >>> filtering cellSeg geometries to match 6272 cells with counts > 0
#> >>> filtering nucSeg geometries to match 6158 cells with counts > 0
#> >>> Reading transcript coordinates
#> >>> Total of 116 features/genes with no transcript detected or `min_phred` < 20 are removed from SFE object
#> >>> To keep all features -> set `min_phred = NULL`
rowGeometries(sfe)
#> List of length 1
#> names(1): txSpots
rowGeometryNames(sfe)
#> [1] "txSpots"
tx <- rowGeometry(sfe, "txSpots")
txSpots(sfe)
#> Simple feature collection with 398 features and 2 fields
#> Geometry type: MULTIPOINT
#> Dimension:     XYZ
#> Bounding box:  xmin: 0.001159668 ymin: -1008.098 xmax: 1026.159 ymax: -4.473946
#> z_range:       zmin: 13.56251 zmax: 27.21318
#> CRS:           NA
#> First 10 features:
#>            gene codeword_index                       geometry
#> ABCC11   ABCC11             87 MULTIPOINT Z ((111.2996 -47...
#> ACE2       ACE2             31 MULTIPOINT Z ((204.9438 -21...
#> ACKR1     ACKR1            349 MULTIPOINT Z ((8.244019 -88...
#> ACTA2     ACTA2            342 MULTIPOINT Z ((4.195129 -40...
#> ACTG2     ACTG2            231 MULTIPOINT Z ((6.257874 -30...
#> ADAM28   ADAM28            119 MULTIPOINT Z ((0.2313232 -4...
#> ADAMTS1 ADAMTS1            242 MULTIPOINT Z ((5.13385 -321...
#> ADGRE1   ADGRE1             24 MULTIPOINT Z ((119.8132 -7....
#> ADGRL4   ADGRL4            132 MULTIPOINT Z ((4.375183 -32...
#> ADH1C     ADH1C             92 MULTIPOINT Z ((1.244507 -59...
unlink(dir_use, recursive = TRUE)
```
