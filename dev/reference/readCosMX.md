# Read CosMX data into SFE

This function reads the standard CosMX output into an SFE object, as in
"Basic Data Files" on the Nanostring website. For new version of CosMX,
these files are the flat files in the AtoMX output.

## Usage

``` r
readCosMX(
  data_dir,
  z = "all",
  sample_id = "sample01",
  min_area = NULL,
  add_molecules = FALSE,
  split_cell_comps = FALSE,
  BPPARAM = SerialParam(),
  file_out = file.path(data_dir, "tx_spots.parquet"),
  z_option = c("3d", "split")
)
```

## Arguments

- data_dir:

  Top level output directory.

- z:

  Integer z index or "all" to indicate which z-planes to read for the
  transcript spots.

- sample_id:

  A `character` sample identifier, which matches the `sample_id` in
  [`imgData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html).
  The `sample_id` will also be stored in a new column in
  [`colData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-colData.html),
  if not already present. Default = `sample01`.

- min_area:

  Minimum cell area in square microns or pixel units (eg for CosMX).
  Anything smaller will be considered artifact or debris and removed.
  Default to \`NULL\`, ie no filtering of polygons.

- add_molecules:

  Logical, whether to add transcripts coordinates to an object.

- split_cell_comps:

  Logical, whether to split transcript spot geometries by cell
  compartment. Only relevant when \`add_molecules = TRUE\`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object specifying parallel processing backend and number of threads to
  use for parallelizable tasks:

  1.  To load cell segmentation from HDF5 files from different fields of
      view (FOVs) with multiple cores. A progress bar can be configured
      in the
      [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
      object. When there are numerous FOVs, reading in the geometries
      can be time consuming, so we recommend using a server and larger
      number of threads. This argument is not used if
      `use_cellpose = TRUE` and the parquet file is present.

  2.  To get the largest piece and see if it's larger than `min_area`
      when there are multiple pieces in the cell segmentation for one
      cell.

- file_out:

  Name of file to save the geometry or raster to disk. Especially when
  the geometries are so large that it's unwieldy to load everything into
  memory. If this file (or directory for multiple files) already exists,
  then the existing file(s) will be read, skipping the processing. When
  writing the file, extensions supplied are ignored and extensions are
  determined based on \`dest\`.

- z_option:

  What to do with z coordinates. "3d" is to construct 3D geometries.
  "split" is to create a separate 2D geometry for each z-plane so
  geometric operations are fully supported but some data wrangling is
  required to perform 3D analyses. When the z coordinates are not
  integers, 3D geometries will always be constructed since there are no
  z-planes to speak of. This argument does not apply when
  \`spatialCoordsNames\` has length 2.

## Value

An SFE object. Cell polygons are written to
\`cell_boundaries_sf.parquet\` in \`data_dir\`. If reading transcript
spots (\`add_molecules = TRUE\`), then the reformatted transcript spots
are saved to file specified in the \`file_out\` argument, which is by
default \`tx_spots.parquet\` in the same directory as the rest of the
data.

## Examples

``` r
fp <- tempfile()
dir_use <- SFEData::CosMXOutput(file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> The downloaded files are in /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce86d1d730b/cosmx 
sfe <- readCosMX(dir_use, z = "all", add_molecules = TRUE)
#> >>> Constructing cell polygons
#> >>> Checking polygon validity
#> >>> Reading transcript coordinates
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk
# Clean up
unlink(dir_use, recursive = TRUE)
```
