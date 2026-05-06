# Read 10X Visium data as SpatialFeatureExperiment

Read Space Ranger output from Visium v1 (not HD) as a
SpatialFeatureExperiment object, where spots are represented with
polygons in the colGeometry called "spotPoly". Other geometries can be
added later after the dataset is read. If `data = "filtered"`, then
spatial neighborhood graphs of the spots are also computed and stored in
the colGraph called "visium" in all samples for downstream spatial
analyses.

## Usage

``` r
read10xVisiumSFE(
  dirs,
  sample_id = paste0("sample", sprintf("%02d", seq_along(dirs))),
  type = c("HDF5", "sparse"),
  data = c("filtered", "raw"),
  images = c("lowres", "hires"),
  unit = c("full_res_image_pixel", "micron"),
  style = "W",
  zero.policy = NULL,
  row.names = c("id", "symbol"),
  flip = c("geometry", "image", "none"),
  read_spatial_enrichment = TRUE
)
```

## Arguments

- dirs:

  Directory for each sample that contains the `spatial` and
  `raw/filtered_featues_bc_matrix` directories.

- sample_id:

  Which sample(s) in the SFE object to use for the graph. Can also be
  "all", which means this function will compute the graph for all
  samples independently.

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

- read_spatial_enrichment:

  Logical, whether to read the \`spatial_enrichment.csv\` file from
  Visium output if the file is present and add its contents to
  \`rowData\`.

## Value

A SpatialFeatureExperiment object. The images might need to be manually
transposed and/or mirrored to match the spots in this version of this
package.

## Note

It is assumed that the images have not been cropped. Otherwise the
images might not align with the spots.

## Examples

``` r
dir <- system.file("extdata", package = "SpatialFeatureExperiment")

sample_ids <- c("sample01", "sample02")
samples <- file.path(dir, sample_ids, "outs")

list.files(samples[1])
#> [1] "filtered_feature_bc_matrix" "spatial"                   
list.files(file.path(samples[1], "spatial"))
#> [1] "aligned_fiducials.jpg"              "barcode_fluorescence_intensity.csv"
#> [3] "detected_tissue_image.jpg"          "scalefactors_json.json"            
#> [5] "spatial_enrichment.csv"             "tissue_hires_image.png"            
#> [7] "tissue_lowres_image.png"            "tissue_positions.csv"              
(sfe <- read10xVisiumSFE(dirs = samples, sample_id = sample_ids,
    type = "sparse", data = "filtered"
))
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample01
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample02
#> class: SpatialFeatureExperiment 
#> dim: 5 25 
#> metadata(0):
#> assays(1): counts
#> rownames(5): ENSG00000014257 ENSG00000142515 ENSG00000263639
#>   ENSG00000163810 ENSG00000149591
#> rowData names(14): symbol Feature.Type ...
#>   Median.Normalized.Average.Counts_sample02
#>   Barcodes.Detected.per.Feature_sample02
#> colnames(25): GTGGCGTGCACCAGAG-1 GGTCCCATAACATAGA-1 ...
#>   TGCAATTTGGGCACGG-1 ATGCCAATCGCTCTGC-1
#> colData names(10): in_tissue array_row ... channel3_mean channel3_stdev
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
#> 
#> unit: full_res_image_pixel
#> Geometries:
#> colGeometries: spotPoly (POLYGON) 
#> 
#> Graphs:
#> sample01: col: visium
#> sample02: col: visium
```
