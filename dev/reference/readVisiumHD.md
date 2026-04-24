# Read Visium HD data

This function reads Visium HD Space Ranger output into R.

## Usage

``` r
readVisiumHD(
  data_dir,
  bin_size = 8L,
  sample_id = "sample01",
  use_cellseg = TRUE,
  type = c("HDF5", "sparse"),
  data = c("filtered", "raw"),
  images = c("lowres", "hires"),
  unit = c("full_res_image_pixel", "micron"),
  style = "W",
  zero.policy = TRUE,
  row.names = c("id", "symbol"),
  flip = c("geometry", "image"),
  add_graph = FALSE
)
```

## Arguments

- data_dir:

  Directory with Visium HD output

- bin_size:

  One or more resolutions to load, must be 2, 8, or 16. Can be either
  integer or character.

- sample_id:

  Which sample(s) in the SFE object to use for the graph. Can also be
  "all", which means this function will compute the graph for all
  samples independently.

- use_cellseg:

  Logical, whether to use cell segmentation if present. Only applicable
  when the `segmented_outputs` directory is present. If TRUE and the
  directory is present, then `bin_size` will be ignored and cell
  segmentations will be used instead.

- type:

  Either "HDF5", and the matrix will be represented as `TENxMatrix`, or
  "sparse", and the matrix will be read as a `dgCMatrix`.

- data:

  character string specifying whether to read in filtered (spots mapped
  to tissue) or raw data (all spots).

- images:

  character vector specifying which images to include. Valid values are
  `"lowres", "hires", "fullres", "detected", "aligned"`

- unit:

  Whether to use pixels in full resolution image or microns as the unit.
  If using microns, then spacing between spots in pixels will be used to
  convert the coordinates into microns, as the spacing is known to be
  100 microns. This is used to plot scale bar.

- style:

  `style` can take values “W”, “B”, “C”, “U”, “minmax” and “S”

- zero.policy:

  default NULL, use global option value; if FALSE stop with error for
  any empty neighbour sets, if TRUE permit the weights list to be formed
  with zero-length weights vectors

- row.names:

  String specifying whether to use Ensembl IDs ("ID") or gene symbols
  ("Symbol") as row names. For symbols, the Ensembl ID will be appended
  to disambiguate rows where the same symbol corresponds to multiple
  Ensembl IDs.

- flip:

  Whether to flip the geometries or the images, because in `sf` and
  `terra`, the geometries use the Cartesian coordinates greater y
  coordinates going up, while in images, greater y values go down.
  Originally the Visium spots are in pixels in full res image. Either
  the image or the geometry needs to be flipped for them match in the
  Cartesian coordinate system.

- add_graph:

  Logical, if to add spatial neighborhood graph for spots and only if
  `c(data = "filtered")`. Default is `FALSE`. This is optional because
  for larger datasets, the graph can take a while to compute.

## Value

An SFE object if \`length(bin_size) == 1L\`, otherwise a list of SFE
objects each element of which is for one bin size. They're not
concatenated since it might not make sense to perform joint analyses on
the different resolutions that benefit from having them in the same SFE
object, unlike different biological replica. Here unlike in
[`read10xVisiumSFE`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/read10xVisiumSFE.md),
the centroids geometry is also added because it will greatly facilitate
plotting when there are many spots when not zooming in. See the
`scattermore` argument in
[`plotSpatialFeature`](https://rdrr.io/pkg/Voyager/man/plotSpatialFeature.html).

## Examples

``` r
library(OSTA.data)
id <- "VisiumHD_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
#> Requesting folder 'segmented_outputs' from OSF
#> Downloaded 6 file(s) from OSF folder 'segmented_outputs'
#> Requesting folder 'binned_outputs' from OSF
#> Downloaded 11 file(s) from OSF folder 'binned_outputs'
dir.create(dir <- tempfile())
unzip(pa, exdir=dir)

# Use binned output
sfe1 <- readVisiumHD(dir, bin_size = 16, use_cellseg = FALSE)
#> >>> 10X VisiumHD data will be loaded: square_016um
#> Skipping missing images
#>   /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce814eeeade/binned_outputs/square_016um/spatial/tissue_hires_image.png

# Use cell segmentations
sfe <- readVisiumHD(dir, use_cellseg = TRUE, unit = "micron")
#> >>> 10X VisiumHD data will be loaded: segmented_outputs
#> >>> Converting pixels to microns

unlink(dir, recursive = TRUE)
```
