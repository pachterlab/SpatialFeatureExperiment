# Get relevant fields and file paths for transcript spots

Get column names for x, y, and z coordinates, gene IDs, and cell IDs
from the transcript file and get file paths for transcript spot
coordinates given technology.

## Usage

``` r
getTechTxFields(tech, data_dir = NULL)
```

## Arguments

- tech:

  Name of the commercial technology, must be one of Vizgen, Xenium, and
  CosMX.

- data_dir:

  Top level directory of the output.

## Value

A named list with elements:

- `spatialCoordsNames`:

  A character vector for column names for the xyz coordinates of the
  transcript spots.

- `gene_col`:

  Column name for gene IDs.

- `cell_col`:

  Column name for cell IDs.

- `fn`:

  File path of the transcript spot file.
