# Aggregate transcript spots from file

This function reads the transcript spot file from the standard output of
the commercial technologies (not GeoParquet) for spatial aggregation
where the spots are assigned to polygons such as cells or spatial bins.
Presets for Xenium, MERFISH, and CosMX are available. For Vizgen and
Xenium, the images can be added when `add_images = TRUE`.

## Usage

``` r
aggregateTx(
  file,
  df = NULL,
  by = NULL,
  sample_id = "sample01",
  spatialCoordsNames = c("X", "Y", "Z"),
  gene_col = "gene",
  phred_col = "qv",
  min_phred = 20,
  flip_geometry = FALSE,
  cellsize = NULL,
  square = TRUE,
  flat_topped = FALSE,
  new_geometry_name = "bins",
  unit = "micron",
  sparse = FALSE,
  BPPARAM = SerialParam(),
  save_memory = FALSE,
  .orig_nrows = NULL
)

aggregateTxTech(
  data_dir,
  df = NULL,
  by = NULL,
  tech = c("Vizgen", "Xenium", "CosMX"),
  sample_id = "sample01",
  image = NULL,
  min_phred = 20,
  flip = c("geometry", "image", "none"),
  max_flip = "50 MB",
  cellsize = NULL,
  square = TRUE,
  flat_topped = FALSE,
  new_geometry_name = "bins",
  sparse = FALSE,
  BPPARAM = SerialParam(),
  save_memory = FALSE
)
```

## Arguments

- file:

  File with the transcript spot coordinates. Should be one row per spot
  when read into R and should have columns for coordinates on each axis,
  gene the transcript is assigned to, and optionally cell the transcript
  is assigned to. Must be csv, tsv, or parquet.

- df:

  If the file is already loaded into memory, a data frame (sf) with
  columns for the x, y, and optionally z coordinates and gene assignment
  of each transcript spot. If specified, then argument `file` will be
  ignored.

- by:

  A `sfc` or `sf` object for spatial aggregation.

- sample_id:

  Which sample in the SFE object the transcript spots should be added
  to.

- spatialCoordsNames:

  Column names for the x, y, and optionally z coordinates of the spots.
  The defaults are for Vizgen.

- gene_col:

  Column name for genes.

- phred_col:

  Column name for Phred scores of the spots.

- min_phred:

  Minimum Phred score to keep spot. By default 20, the conventional
  threshold indicating "acceptable", meaning that there's 1 chance that
  the spot was decoded in error.

- flip_geometry:

  Logical, whether to flip the transcript spot geometries to match the
  images if added later.

- cellsize:

  numeric of length 1 or 2 with target cellsize: for square or
  rectangular cells the width and height, for hexagonal cells the
  distance between opposite edges (edge length is cellsize/sqrt(3)). A
  length units object can be passed, or an area unit object with area
  size of the square or hexagonal cell.

- square:

  logical; if `FALSE`, create hexagonal grid

- flat_topped:

  logical; if `TRUE` generate flat topped hexagons, else generate pointy
  topped

- new_geometry_name:

  Name to give to the new `colGeometry` in the output. Defaults to
  "bins".

- unit:

  Unit the coordinates are in, either microns or pixels in full
  resolution image.

- sparse:

  Logical, whether the gene count matrix from aggregating transcript
  spots should be sparse. When the bins are large, the matrix will not
  be very sparse so using sparse matrix will not save memory, but when
  the bins are small, sparsity is worth it.

- BPPARAM:

  bpparam object to specify parallel computing over genes. If a lot of
  memory is used, then stick to \`SerialParam()\`. Using multicore in
  the R level when \`save_memory = TRUE\` works and results into faster
  computation.

- save_memory:

  Logical, if TRUE, then the transcript spots will not all be loaded
  into memory.
  [`open_dataset`](https://arrow.apache.org/docs/r/reference/open_dataset.html)
  is used to open a link to the data and then transcript spots of one
  gene is loaded into memory at a time.

- .orig_nrows:

  Only used internally in the SFE method of `aggregate`

- data_dir:

  Top level output directory.

- tech:

  Which technology whose output to read, must be one of "Vizgen",
  "Xenium", or "CosMX" though more technologies may be added later.

- image:

  String, which image(s) to add to the output SFE object. Not applicable
  to CosMX. See
  [`readVizgen`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readVizgen.md)
  and
  [`readXenium`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readXenium.md)
  for options and multiple images can be specified. If `NULL`, then the
  default from the read function for the technology will be used.

- flip:

  Logical, whether to flip the geometry to match image. Here the y
  coordinates are simply set to -y, so the original bounding box is not
  preserved. This is consistent with `readVizgen` and `readXenium`.

- max_flip:

  Maximum size of the image allowed to flip the image. Because the image
  will be loaded into memory to be flipped. If the image is larger than
  this size then the coordinates will be flipped instead.

## Value

A SFE object with count matrix for number of spots of each gene in each
geometry. Geometries with no spot are removed.

## Note

The resulting SFE object often includes geometries (e.g. grid cells)
outside tissue, because there can be transcript spots detected outside
the tissue. Also, bins at the edge of the tissue that don't fully
overlap with the tissue will have lower transcript counts; this may have
implications to downstream spatial analyses.
