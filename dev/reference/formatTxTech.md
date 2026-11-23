# Read and process transcript spots for specific commercial technologies

To preset parameters such as `spatialCoordsNames`, `gene_col`,
`cell_col`, and `phred_col` that are standard for the output of the
technology.

## Usage

``` r
formatTxTech(
  data_dir,
  tech = c("Vizgen", "Xenium", "CosMX"),
  dest = c("rowGeometry", "colGeometry"),
  z = deprecated(),
  min_phred = 20,
  split_cell_comps = deprecated(),
  z_option = deprecated(),
  flip = FALSE,
  file_out = NULL,
  BPPARAM = SerialParam(),
  return = TRUE
)

addTxTech(
  sfe,
  data_dir,
  sample_id = 1L,
  tech = c("Vizgen", "Xenium", "CosMX"),
  z = deprecated(),
  min_phred = 20,
  split_cell_comps = deprecated(),
  z_option = deprecated(),
  flip = FALSE,
  file_out = NULL,
  BPPARAM = SerialParam()
)
```

## Arguments

- data_dir:

  Top level output directory.

- tech:

  Which technology whose output to read, must be one of "Vizgen",
  "Xenium", or "CosMX" though more technologies may be added later.

- dest:

  This argument and the "colGeometry" option is deprecated. Only
  "rowGeometry" will be allowed.

- z:

  Which z-planes to read. Always "all" for Xenium where the z
  coordinates are not discrete.

- min_phred:

  Minimum Phred score to keep spot. By default 20, the conventional
  threshold indicating "acceptable", meaning that there's 1% chance that
  the spot was decoded in error.

- split_cell_comps:

  Deprecated.

- z_option:

  Deprecated. Will always be "3d" to construct 3D point geometries if
  multiple z-planes are present.

- flip:

  Logical, whether to flip the geometry to match image. Here the y
  coordinates are simply set to -y, so the original bounding box is not
  preserved. This is consistent with `readVizgen` and `readXenium`.

- file_out:

  Name of file to save the geometry or raster to disk. Especially when
  the geometries are so large that it's unwieldy to load everything into
  memory. If this file (or directory for multiple files) already exists,
  then the existing file(s) will be read, skipping the processing. When
  writing the file, extensions supplied are ignored and extensions are
  determined based on `dest`.

- BPPARAM:

  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object to specify multithreading to convert raw char in some parquet
  files to R objects. Not used otherwise.

- return:

  Logical, whether to return the geometries in memory. This does not
  depend on whether the geometries are written to file. Always `FALSE`
  when `dest = "colGeometry"`.

- sfe:

  A `SpatialFeatureExperiment` object.

- sample_id:

  Which sample in the SFE object the transcript spots should be added
  to.

## Value

The `sf` data frame, or path to file where geometries are written if
`return = FALSE`.

## Examples

``` r
library(SFEData)
fp <- tempfile()
dir_use <- XeniumOutput("v2", file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> The downloaded files are in /private/var/folders/q2/p1t6rgxj5fj32g1dc7dmd49r0000gn/T/Rtmp4UEbWg/file1808404a2270/xenium2 
fn_tx <- formatTxTech(dir_use, tech = "Xenium", flip = TRUE, return = FALSE,
                      file_out = file.path(dir_use, "tx_spots.parquet"))
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk
```
