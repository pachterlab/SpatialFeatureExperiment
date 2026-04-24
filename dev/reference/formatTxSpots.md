# Read and process transcript spots geometry for SFE

The function \`formatTxSpots\` reads the transcript spot coordinates of
smFISH-based data and formats the data. The data is not added to an SFE
object. If the file specified in \`file_out\` already exists, then this
file will be read instead of the original file in the \`file\` argument,
so the processing is not run multiple times. The function \`addTxSpots\`
adds the data read and processed in \`formatTxSpots\` to the SFE object,
and reads all transcript spot data. To only read a subset of transcript
spot data, first use \`formatTxSpots\` to write the re-formatted data to
disk. Then read the specific subset and add them separately to the SFE
object with the setter functions.

## Usage

``` r
formatTxSpots(
  file,
  dest = c("rowGeometry", "colGeometry"),
  spatialCoordsNames = c("global_x", "global_y", "global_z"),
  gene_col = "gene",
  cell_col = "cell_id",
  z = "all",
  phred_col = "qv",
  min_phred = 20,
  split_col = NULL,
  not_in_cell_id = c("-1", "UNASSIGNED"),
  z_option = c("3d", "split"),
  flip = FALSE,
  file_out = NULL,
  BPPARAM = SerialParam(),
  return = TRUE,
  save_memory = FALSE,
  partition = FALSE
)

addTxSpots(
  sfe,
  file,
  sample_id = 1L,
  spatialCoordsNames = c("global_x", "global_y", "global_z"),
  gene_col = "gene",
  z = "all",
  phred_col = "qv",
  min_phred = 20,
  split_col = NULL,
  z_option = c("3d", "split"),
  flip = FALSE,
  file_out = NULL,
  BPPARAM = SerialParam()
)
```

## Arguments

- file:

  File with the transcript spot coordinates. Should be one row per spot
  when read into R and should have columns for coordinates on each axis,
  gene the transcript is assigned to, and optionally cell the transcript
  is assigned to. Must be csv, tsv, or parquet.

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

- spatialCoordsNames:

  Column names for the x, y, and optionally z coordinates of the spots.
  The defaults are for Vizgen.

- gene_col:

  Column name for genes.

- cell_col:

  Column name for cell IDs, ignored if \`dest = "rowGeometry"\`. Can
  have length \> 1 when multiple columns are needed to uniquely identify
  cells, in which case the contents of the columns will be concatenated,
  such as in CosMX data where cell ID is only unique within the same
  FOV. Default "cell_id" is for Vizgen MERFISH. Should be \`c("cell_ID",
  "fov")\` for CosMX.

- z:

  Index of z plane to read. Can be "all" to read all z-planes into
  MULTIPOINT geometries with XYZ coordinates. If z values are not
  integer, then spots with all z values will be read.

- phred_col:

  Column name for Phred scores of the spots.

- min_phred:

  Minimum Phred score to keep spot. By default 20, the conventional
  threshold indicating "acceptable", meaning that there's 1 chance that
  the spot was decoded in error.

- split_col:

  Categorical column to split the geometries, such as cell compartment
  the spots are assigned to as in the "CellComp" column in CosMX output.

- not_in_cell_id:

  Value of cell ID indicating that the spot is not assigned to any cell,
  such as "-1" in Vizgen MERFISH and "0" in CosMX. When there're
  multiple columns for \`cell_col\`, the first column is used to
  identify spots that are not in cells.

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

- save_memory:

  Don't load transcript spots of all genes at once. Not fully
  implemented here.

- partition:

  Whether to partition the output by gene.

- sfe:

  A \`SpatialFeatureExperiment\` object.

- sample_id:

  Which sample in the SFE object the transcript spots should be added
  to.

## Value

A sf data frame for vector geometries if \`file_out\` is not set.
\`SpatRaster\` for raster. If there are multiple files written, such as
when splitting by cell compartment or when \`dest = "colGeometry"\`,
then a directory with the same name as \`file_out\` will be created (but
without the extension) and the files are written to that directory with
informative names. \`parquet\` files that can be read with \`st_read\`
is written for vector geometries. When \`return = FALSE\`, the file name
or directory (when there're multiple files) is returned.

The \`sf\` data frame, or path to file where geometries are written if
\`return = FALSE\`.

## Note

When \`dest = "colGeometry"\`, the geometries are always written to disk
and not returned in memory, because this is essentially the gene count
matrix, which is sparse. This kind of reformatting is implemented so
users can read in MULTIPOINT geometries with transcript spots for each
gene assigned to each cell for spatial point process analyses, where not
all genes are loaded at once.

## Examples

``` r
# Default arguments are for MERFISH
fp <- tempfile()
dir_use <- SFEData::VizgenOutput(file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> The downloaded files are in /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce8772603dd/vizgen_cellbound 
g <- formatTxSpots(file.path(dir_use, "detected_transcripts.csv"))
#> >>> Converting transcript spots to geometry
unlink(dir_use, recursive = TRUE)

# For CosMX, note the colnames, also dest = "colGeometry"
# Results are written to the tx_spots directory
dir_use <- SFEData::CosMXOutput(file_path = fp)
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> The downloaded files are in /private/var/folders/tb/y368xp_x10s3ty1b_mtl5mxr0000gn/T/RtmpNrtF20/file1ce8772603dd/cosmx 
cg <- formatTxSpots(file.path(dir_use, "Run5642_S3_Quarter_tx_file.csv"),
dest = "colGeometry", z = "all",
cell_col = c("cell_ID", "fov"),
gene_col = "target", not_in_cell_id = "0",
spatialCoordsNames = c("x_global_px", "y_global_px", "z"),
file_out = file.path(dir_use, "tx_spots"))
#> >>> Converting transcript spots to geometry
#> >>> Writing reformatted transcript spots to disk
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |                                                                      |   1%
  |                                                                            
  |=                                                                     |   1%
  |                                                                            
  |=                                                                     |   2%
  |                                                                            
  |==                                                                    |   2%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |===                                                                   |   4%
  |                                                                            
  |===                                                                   |   5%
  |                                                                            
  |====                                                                  |   5%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |=====                                                                 |   7%
  |                                                                            
  |=====                                                                 |   8%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |======                                                                |   9%
  |                                                                            
  |=======                                                               |   9%
  |                                                                            
  |=======                                                               |  10%
  |                                                                            
  |=======                                                               |  11%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |========                                                              |  12%
  |                                                                            
  |=========                                                             |  12%
  |                                                                            
  |=========                                                             |  13%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |==========                                                            |  15%
  |                                                                            
  |===========                                                           |  15%
  |                                                                            
  |===========                                                           |  16%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |============                                                          |  18%
  |                                                                            
  |=============                                                         |  18%
  |                                                                            
  |=============                                                         |  19%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |==============                                                        |  20%
  |                                                                            
  |==============                                                        |  21%
  |                                                                            
  |===============                                                       |  21%
  |                                                                            
  |===============                                                       |  22%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |================                                                      |  23%
  |                                                                            
  |=================                                                     |  24%
  |                                                                            
  |=================                                                     |  25%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |==================                                                    |  26%
  |                                                                            
  |===================                                                   |  27%
  |                                                                            
  |===================                                                   |  28%
  |                                                                            
  |====================                                                  |  28%
  |                                                                            
  |====================                                                  |  29%
  |                                                                            
  |=====================                                                 |  29%
  |                                                                            
  |=====================                                                 |  30%
  |                                                                            
  |=====================                                                 |  31%
  |                                                                            
  |======================                                                |  31%
  |                                                                            
  |======================                                                |  32%
  |                                                                            
  |=======================                                               |  32%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |========================                                              |  34%
  |                                                                            
  |========================                                              |  35%
  |                                                                            
  |=========================                                             |  35%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |==========================                                            |  37%
  |                                                                            
  |==========================                                            |  38%
  |                                                                            
  |===========================                                           |  38%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |============================                                          |  39%
  |                                                                            
  |============================                                          |  40%
  |                                                                            
  |============================                                          |  41%
  |                                                                            
  |=============================                                         |  41%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |==============================                                        |  42%
  |                                                                            
  |==============================                                        |  43%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |===============================                                       |  45%
  |                                                                            
  |================================                                      |  45%
  |                                                                            
  |================================                                      |  46%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |=================================                                     |  48%
  |                                                                            
  |==================================                                    |  48%
  |                                                                            
  |==================================                                    |  49%
  |                                                                            
  |===================================                                   |  49%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |===================================                                   |  51%
  |                                                                            
  |====================================                                  |  51%
  |                                                                            
  |====================================                                  |  52%
  |                                                                            
  |=====================================                                 |  52%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |======================================                                |  54%
  |                                                                            
  |======================================                                |  55%
  |                                                                            
  |=======================================                               |  55%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |========================================                              |  57%
  |                                                                            
  |========================================                              |  58%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |=========================================                             |  59%
  |                                                                            
  |==========================================                            |  59%
  |                                                                            
  |==========================================                            |  60%
  |                                                                            
  |==========================================                            |  61%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |===========================================                           |  62%
  |                                                                            
  |============================================                          |  62%
  |                                                                            
  |============================================                          |  63%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |=============================================                         |  65%
  |                                                                            
  |==============================================                        |  65%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |===============================================                       |  68%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |================================================                      |  69%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |=================================================                     |  70%
  |                                                                            
  |=================================================                     |  71%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |==================================================                    |  72%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |===================================================                   |  73%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |=====================================================                 |  75%
  |                                                                            
  |=====================================================                 |  76%
  |                                                                            
  |======================================================                |  77%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |=======================================================               |  78%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |========================================================              |  79%
  |                                                                            
  |========================================================              |  80%
  |                                                                            
  |========================================================              |  81%
  |                                                                            
  |=========================================================             |  81%
  |                                                                            
  |=========================================================             |  82%
  |                                                                            
  |==========================================================            |  82%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |===========================================================           |  84%
  |                                                                            
  |===========================================================           |  85%
  |                                                                            
  |============================================================          |  85%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |=============================================================         |  87%
  |                                                                            
  |=============================================================         |  88%
  |                                                                            
  |==============================================================        |  88%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |===============================================================       |  89%
  |                                                                            
  |===============================================================       |  90%
  |                                                                            
  |===============================================================       |  91%
  |                                                                            
  |================================================================      |  91%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |=================================================================     |  92%
  |                                                                            
  |=================================================================     |  93%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |==================================================================    |  95%
  |                                                                            
  |===================================================================   |  95%
  |                                                                            
  |===================================================================   |  96%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |====================================================================  |  98%
  |                                                                            
  |===================================================================== |  98%
  |                                                                            
  |===================================================================== |  99%
  |                                                                            
  |======================================================================|  99%
  |                                                                            
  |======================================================================| 100%
#> 
# Cleanup
unlink(dir_use, recursive = TRUE)
```
