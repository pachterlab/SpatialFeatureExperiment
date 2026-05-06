# Read 10X Xenium output as SpatialFeatureExperiment

This function reads the standard 10X Xenium output into an SFE object.

## Usage

``` r
readXenium(
  data_dir,
  sample_id = "sample01",
  min_area = NULL,
  image = c("morphology_focus", "morphology_mip"),
  segmentations = c("cell", "nucleus"),
  row.names = c("id", "symbol"),
  flip = c("geometry", "image", "none"),
  max_flip = "50 MB",
  filter_counts = FALSE,
  add_molecules = FALSE,
  min_phred = 20,
  BPPARAM = SerialParam(),
  file_out = file.path(data_dir, "tx_spots.parquet")
)
```

## Arguments

- data_dir:

  Top level output directory.

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

  Which image(s) to load, can be "morphology_mip", "morphology_focus" or
  both. Note that in Xenium Onboarding Analysis (XOA) v2, there is no
  longer "morphology_mip" and "morphology_focus" is a directory with 4
  images corresponding to 4 channels: DAPI, "Cadherin", 18S, and
  Vimentin. So this argument is ignored for XOA v2.

- segmentations:

  Which segmentation outputs to read, can be "cell", "nucleus", or both.

- row.names:

  String specifying whether to use Ensembl IDs ("id") or gene symbols
  ("symbol") as row names. If using symbols, the Ensembl ID will be
  appended to disambiguate in case the same symbol corresponds to
  multiple Ensembl IDs. Always "symbol" if \`add_molecules = TRUE\`
  because only gene symbols are used in the transcript spot files.

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

- min_phred:

  Minimum Phred score to keep spot. By default 20, the conventional
  threshold indicating "acceptable", meaning that there's 1 chance that
  the spot was decoded in error.

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

## Value

An SFE object. If reading segmentations, the cell or nuclei segmentation
will be saved to \`cell_boundaries_sf.parquet\` and
\`nucleus_boundaries_sf.parquet\` respectively in \`data.dir\` so next
time the boundaries can be read much more quickly. If reading transcript
spots (\`add_molecules = TRUE\`), then the reformatted transcript spots
are saved to file specified in the \`file_out\` argument, which is by
default \`tx_spots.parquet\` in the same directory as the rest of the
data. If images are present, then the images will be of the
`BioFormatsImage` class and not loaded into memory until necessary in
later operations.

## Note

Sometimes when reading images, you will see this error the first time:
'java.lang.NullPointerException: Cannot invoke
"loci.formats.DimensionSwapper.setMetadataFiltered(boolean)" because
"RBioFormats.reader" is null'. See this issue
https://github.com/aoles/RBioFormats/issues/42 Rerun the code and it
should work the second time.

## Examples

``` r
library(SFEData)
library(RBioFormats)
fp <- tempfile()
dir_use <- XeniumOutput("v2", file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
#> The downloaded files are in /tmp/RtmphTOllm/fileea9d2e2c821a/xenium2 
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
unlink(dir_use, recursive = TRUE)
```
