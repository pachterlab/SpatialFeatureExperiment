# Save SpatialFeatureExperiment as RDS file

Saving SFE objects as RDS files is complicated by the `SpatRaster` class
of the images. If present, the images need to be wrapped with the
[`wrap`](https://rspatial.github.io/terra/reference/wrap.html) function
in `terra` before serializing the SFE object. Otherwise the images will
be invalid pointers when the RDS is reloaded. If the image does not fit
in memory and its file source is unknown, then it will be written to a
temporary file, which is reloaded when the RDS file is loaded. When an
SFE object with images is read from an RDS file, the images will not be
unwrapped until necessary.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
saveRDS(
  object,
  file = "",
  ascii = FALSE,
  version = NULL,
  compress = TRUE,
  refhook = NULL
)
```

## Arguments

- object:

  A `SpatialFeatureExperiment` object.

- file:

  a [connection](https://rdrr.io/r/base/connections.html) or the path
  name of the file where the R object is saved to or read from.

- ascii:

  a logical. If `TRUE` or `NA`, an ASCII representation is written;
  otherwise (default), a binary one is used. See the comments in the
  help for [`save`](https://rdrr.io/r/base/save.html).

- version:

  the workspace format version to use. `NULL` specifies the current
  default version (3). The only other supported value is 2, the default
  from R 1.4.0 to R 3.5.0.

- compress:

  a logical specifying whether saving to a named file is to use `"gzip"`
  compression, or one of `"gzip"`, `"bzip2"`, `"xz"` or `"zstd"` to
  indicate the type of compression to be used. Ignored if `file` is a
  connection.

- refhook:

  a hook function for handling reference objects. See
  [`unserialize`](https://rdrr.io/r/base/serialize.html).

## Value

Invisibly `NULL`.

## Examples

``` r
outdir <- system.file("extdata", package = "SpatialFeatureExperiment")
samples <- file.path(outdir, paste0("sample0", 1:2))
sfe <- read10xVisiumSFE(samples, type = "sparse", data = "filtered")
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample01
#> >>> 10X Visium data will be loaded: outs
#> >>> Adding spatial neighborhood graph to sample02
saveRDS(sfe, "foo.rds")
#> Warning: This object contains out of memory data that may
#>                           break when the RDS file is moved.
# Clean up
file.remove("foo.rds")
#> [1] TRUE
```
