# Read Vizgen MERFISH output as SpatialFeatureExperiment

This function reads the standard Vizgen MERFISH output into an SFE
object. The coordinates are in microns. Cell centroids are read into
[`colGeometry`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
"centroids", and cell segmentations are read into `colGeometry`
"cellSeg". The image(s) (polyT, DAPI, and cell boundaries) are also read
as
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
objects so they are not loaded into memory unless necessary. Because the
image's origin is the top left while the geometry's origin is bottom
left, either the image or the geometry needs to be flipped. Because the
image accompanying MERFISH datasets are usually very large, the
coordinates will be flipped so the flipping operation won't load the
entire image into memory. Large datasets with hundreds of thousands of
cells can take a while to read if reading transcript spots as it takes a
while to convert the spots to MULTIPOINT geometries.

## Usage

``` r
readVizgen(
  data_dir,
  z = "all",
  sample_id = "sample01",
  min_area = NULL,
  image = c("DAPI", "PolyT", "Cellbound"),
  flip = c("geometry", "image", "none"),
  max_flip = "50 MB",
  filter_counts = FALSE,
  add_molecules = FALSE,
  use_bboxes = FALSE,
  use_cellpose = TRUE,
  BPPARAM = SerialParam(),
  file_out = file.path(data_dir, "detected_transcripts.parquet"),
  z_option = c("3d", "split")
)
```

## Arguments

- data_dir:

  Top level output directory.

- z:

  Integer, z index to read, or "all", indicating z-planes of the images
  and transcript spots to read. While cell segmentation seems to have
  multiple z-planes, the segmentation in all z-planes are the same so in
  effect the cell segmentatio is only in 2D.

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

- image:

  Which image(s) to load, can be "DAPI", "PolyT", "Cellbound" or any
  combination of them.

- flip:

  To flip the image, geometry coordinates, or none. Because the image
  has the origin at the top left while the geometry has origin at the
  bottom left, one of them needs to be flipped for them to match. If one
  of them is already flipped, then use "none". The image will not be
  flipped if it's GeoTIFF.

- max_flip:

  Maximum size of the image allowed to flip the image. Because the image
  will be loaded into memory to be flipped. If the image is larger than
  this size then the coordinates will be flipped instead.

- filter_counts:

  Logical, whether to keep cells with counts `> 0`.

- add_molecules:

  Logical, whether to add transcripts coordinates to an object.

- use_bboxes:

  If no segmentation output is present, use `cell_metadata` to make
  bounding boxes instead.

- use_cellpose:

  Whether to read the parquet files from CellPose cell segmentation. If
  `FALSE`, cell segmentation will be read from the HDF5 files. Note that
  reading HDF5 files for numerous FOVs is very slow.

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

A `SpatialFeatureExperiment` object.

## Note

Since the transcript spots file is often very large, we recommend only
using `add_molecules = TRUE` on servers with a lot of memory. If reading
all z-planes, conversion of transcript spot geometry to parquet file
might fail due to arrow data length limit. In a future version, when the
transcript spot geometry is large, it will be written to multiple
separate parquet files which are then concatenated with DuckDB. Also, in
a future version, the transcript spot processing function might be
rewritten in C++ to stream the original CSV file so it's not entirely
loaded into memory.

## Examples

``` r
fp <- tempfile()
dir_use <- SFEData::VizgenOutput(file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> The downloaded files are in /tmp/RtmphTOllm/fileea9d19b7dba1/vizgen_cellbound 
sfe <- readVizgen(dir_use, z = 3L, image = "PolyT",
flip = "geometry")
#> >>> 1 `.parquet` files exist:
#> /tmp/RtmphTOllm/fileea9d19b7dba1/vizgen_cellbound/cell_boundaries.parquet
#> >>> using -> /tmp/RtmphTOllm/fileea9d19b7dba1/vizgen_cellbound/cell_boundaries.parquet
#> >>> Cell segmentations are found in `.parquet` file
#> >>> Checking polygon validity

## Filtering of counts, and addition of molecule coordinates..
sfe <- readVizgen(dir_use, z = 3L, image = "PolyT", filter_counts = TRUE,
add_molecules = TRUE, flip = "geometry")
#> >>> 1 `.parquet` files exist:
#> /tmp/RtmphTOllm/fileea9d19b7dba1/vizgen_cellbound/cell_boundaries.parquet
#> >>> using -> /tmp/RtmphTOllm/fileea9d19b7dba1/vizgen_cellbound/cell_boundaries.parquet
#> >>> Cell segmentations are found in `.parquet` file
#> >>> filtering `cell_metadata` - keep cells with `transcript_count` > 0
#> >>> Checking polygon validity
#> >>> Reading transcript coordinates
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk

unlink(dir_use, recursive = TRUE)
```
