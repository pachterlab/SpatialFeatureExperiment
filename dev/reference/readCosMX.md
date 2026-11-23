# Read CosMX data into SFE

This function reads the standard CosMX output into an SFE object, as in
"Basic Data Files" on the Nanostring website. For new version of CosMX,
these files are the flat files in the AtoMX output.

## Usage

``` r
readCosMX(
  data_dir,
  z = deprecated(),
  sample_id = "sample01",
  min_area = NULL,
  add_molecules = FALSE,
  split_cell_comps = deprecated(),
  BPPARAM = SerialParam(),
  file_out = file.path(data_dir, "tx_spots.parquet"),
  z_option = deprecated()
)
```

## Arguments

- data_dir:

  Top level output directory.

- z:

  Deprecated.

- sample_id:

  A `character` sample identifier, which matches the `sample_id` in
  [`imgData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html).
  The `sample_id` will also be stored in a new column in
  [`colData`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-colData.html),
  if not already present. Default = `sample01`.

- min_area:

  Minimum cell area in square microns or pixel units (eg for CosMX).
  Anything smaller will be considered artifact or debris and removed.
  Default to `NULL`, ie no filtering of polygons.

- add_molecules:

  Logical, whether to add transcripts coordinates to an object.

- split_cell_comps:

  Deprecated.

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
  determined based on `dest`.

- z_option:

  Deprecated. Will always be "3d" to construct 3D point geometries if
  multiple z-planes are present.

## Value

An SFE object. Cell polygons are written to `cell_boundaries_sf.parquet`
in `data_dir`. If reading transcript spots (`add_molecules = TRUE`),
then the reformatted transcript spots are saved to file specified in the
`file_out` argument, which is by default `tx_spots.parquet` in the same
directory as the rest of the data.

## Examples

``` r
fp <- tempfile()
dir_use <- SFEData::CosMXOutput(file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> The downloaded files are in /private/var/folders/q2/p1t6rgxj5fj32g1dc7dmd49r0000gn/T/Rtmp4UEbWg/file18081dc9f2dc/cosmx 
sfe <- readCosMX(dir_use, add_molecules = TRUE)
#> >>> Constructing cell polygons
#> >>> Checking polygon validity
#> >>> Reading transcript coordinates
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk
# Clean up
unlink(dir_use, recursive = TRUE)
```
