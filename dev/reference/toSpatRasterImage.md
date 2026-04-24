# Convert images to SpatRasterImage

The resolution specified from the OME-TIFF file will be read into memory
and written to disk as a GeoTIFF file that has the extent. The output
file will have the same file name as the input file except without the
`ome` in the extension.

## Usage

``` r
# S4 method for class 'ExtImage'
toSpatRasterImage(
  x,
  save_geotiff = TRUE,
  file_out = "img.tiff",
  overwrite = FALSE
)

# S4 method for class 'BioFormatsImage'
toSpatRasterImage(
  x,
  save_geotiff = TRUE,
  resolution = 4L,
  channel = NULL,
  overwrite = FALSE
)
```

## Arguments

- x:

  Either a `BioFormatsImage` or `EBIImage` object.

- save_geotiff:

  Logical, whether to save the image to GeoTIFF file.

- file_out:

  File to save the non-OME TIFF file for `SpatRaster`.

- overwrite:

  Logical, whether to overwrite existing file of the same name.

- resolution:

  Integer, which resolution in the `BioFormatsImage` to read and
  convert. Defaults to 4, which is a lower resolution. Ignored if only 1
  resolution is present.

- channel:

  Integer vector to indicate channel(s) to read. If `NULL`, then all
  channels will be read.

## Value

A `SpatRasterImage` object

## See also

toExtImage
