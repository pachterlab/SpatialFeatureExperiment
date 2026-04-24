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
  z = "all",
  min_phred = 20,
  split_cell_comps = FALSE,
  z_option = c("3d", "split"),
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
  z = "all",
  min_phred = 20,
  split_cell_comps = FALSE,
  z_option = c("3d", "split"),
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

  Where in the SFE object to store the spot geometries. This affects how
  the data is processed. Options:

  rowGeometry

  :   All spots for each gene will be a \`MULTIPOINT\` geometry,
      regardless of whether they are in cells or which cells they are
      assigned to.

  colGeometry

  :   The spots for each gene assigned to a cell of interest will be a
      \`MULTIPOINT\` geometry; since the gene count matrix is sparse,
      the geometries are NOT returned to memory.

- z:

  Which z-planes to read. Always "all" for Xenium where the z
  coordinates are not discrete.

- min_phred:

  Minimum Phred score to keep spot. By default 20, the conventional
  threshold indicating "acceptable", meaning that there's 1 chance that
  the spot was decoded in error.

- split_cell_comps:

  Only relevant to CosMX whose transcript spot file assigns the spots to
  cell components. Setting this argument to `TRUE`

- z_option:

  What to do with z coordinates. "3d" is to construct 3D geometries.
  "split" is to create a separate 2D geometry for each z-plane so
  geometric operations are fully supported but some data wrangling is
  required to perform 3D analyses. When the z coordinates are not
  integers, 3D geometries will always be constructed since there are no
  z-planes to speak of. This argument does not apply when
  \`spatialCoordsNames\` has length 2.

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
  determined based on \`dest\`.

- BPPARAM:

  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object to specify multithreading to convert raw char in some parquet
  files to R objects. Not used otherwise.

- return:

  Logical, whether to return the geometries in memory. This does not
  depend on whether the geometries are written to file. Always \`FALSE\`
  when \`dest = "colGeometry"\`.

- sfe:

  A \`SpatialFeatureExperiment\` object.

- sample_id:

  Which sample in the SFE object the transcript spots should be added
  to.

## Value

The \`sf\` data frame, or path to file where geometries are written if
\`return = FALSE\`.

## Examples

``` r
library(SFEData)
fp <- tempfile()
dir_use <- XeniumOutput("v2", file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> The downloaded files are in /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce8635f7e99/xenium2 
fn_tx <- formatTxTech(dir_use, tech = "Xenium", flip = TRUE, return = FALSE,
                      file_out = file.path(dir_use, "tx_spots.parquet"))
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk
```
