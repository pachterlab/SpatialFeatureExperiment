# Read transcript spots of select genes

I speculate that in practice, the most common use of the transcript
spots is visualization, and only a few genes can be visualized at a time
or the spots will overcrowd. Then it doesn't make sense to load the
transcript spots of all genes into memory as they can take up a lot of
memory. The function `readSelectTx` reads transcript spots of select
genes into R, and the function `addSelectTx` adds them to
`rowGeometries` of the SFE object.

## Usage

``` r
readSelectTx(file, gene_select, z = "all", z_option = c("3d", "split"))

addSelectTx(
  sfe,
  file,
  gene_select,
  sample_id = 1L,
  z = "all",
  z_option = c("3d", "split"),
  swap_rownames = NULL
)
```

## Arguments

- file:

  File path of a GeoParquet file (e.g. already reformatted with the
  [`formatTxSpots`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/formatTxSpots.md)
  or
  [`addTxSpots`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/formatTxSpots.md)
  function, should have already flipped to match image if necessary).

- gene_select:

  Character vector of a subset of genes. If `NULL`, then all genes that
  have transcript spots are added. Only relevant when reading data from
  formatted files on disk. If specified, then `return = TRUE`.

- z:

  Index of z plane to read. Can be "all" to read all z-planes into
  MULTIPOINT geometries with XYZ coordinates. If z values are not
  integer, then spots with all z values will be read.

- z_option:

  What to do with z coordinates. "3d" is to construct 3D geometries.
  "split" is to create a separate 2D geometry for each z-plane so
  geometric operations are fully supported but some data wrangling is
  required to perform 3D analyses. When the z coordinates are not
  integers, 3D geometries will always be constructed since there are no
  z-planes to speak of. This argument does not apply when
  \`spatialCoordsNames\` has length 2.

- sfe:

  A \`SpatialFeatureExperiment\` object.

- sample_id:

  Which sample in the SFE object the transcript spots should be added
  to.

- swap_rownames:

  Name of a column in `rowData(sfe)` to use as gene identifiers in place
  of the actual row names. In some cases this may be needed to match
  each transcript spot MULTIPOINT geometry to rows of `sfe`.

## Value

When there are multipel parquet files to be read, a list of sf data
frames with MULTIPOINT geometry for genes selected. When there is only
one file, then one sf data frame. For `addSelectTx`, an SFE object with
the transcript spots of the selected genes added.

## Note

The GDAL Parquet driver is required for this function, though not for
other functions that work with GeoParquet files. GDAL Parquet driver has
been supported since GDAL 3.5.0, but is not part of the default
installation. The `z` and `z_option` arguments are there since the file
names contain z-plane information when relevant. See the [GDAL
documentation page for the Parquet
driver](https://gdal.org/drivers/vector/parquet.html).

## Examples

``` r
library(SFEData)
if (gdalParquetAvailable()) {
    fp <- tempfile()
    dir_use <- XeniumOutput("v2", file_path = fp)
    fn_tx <- formatTxTech(dir_use, tech = "Xenium", flip = TRUE, return = FALSE,
                          file_out = file.path(dir_use, "tx_spots.parquet"))
    gene_select <- c("ACE2", "BMX")
    df <- readSelectTx(fn_tx, gene_select)
    # RBioFormats null pointer error the first time
    try(sfe <- readXenium(dir_use))
    sfe <- readXenium(dir_use)
    sfe <- addSelectTx(sfe, fn_tx, head(rownames(sfe), 5), swap_rownames = "Symbol")
    unlink(dir_use, recursive = TRUE)
}
```
