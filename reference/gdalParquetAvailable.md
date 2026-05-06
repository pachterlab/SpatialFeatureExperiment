# Check if Parquet GDAL driver is available

The GeoParquet files for geometries are typically written and read with
the `sfarrow` package, but to add only a select few genes to the SFE
object say for visualization purposes, the Parquet GDAL driver is
required in order to use GDAL's SQL to query the GeoParquet file to only
load the few genes requested. The transcript spots from a large dataset
can take up a lot of memory if all loaded.

## Usage

``` r
gdalParquetAvailable()
```

## Value

Logical, indicating whether the Parquet driver is present.

## Details

The Parquet driver has been supported since GDAL 3.5.0. The `arrow` C++
library must be installed in order to make the Parquet driver available.
When arrow is installed, newer versions of GDAL installed from Homebrew
(Mac) should have the Parquet driver. For Linux, the binary from
`apt-get`'s default repo is 3.4.1 (as of April 2024). To use the Parquet
driver, GDAL may need to be installed from source. See script from the
[geospatial
rocker](https://github.com/rocker-org/rocker-versioned2/blob/master/scripts/experimental/install_dev_osgeo.sh).
A Voyager docker container with the Parquet driver will soon be
provided.

## Examples

``` r
gdalParquetAvailable()
#> [1] FALSE
```
