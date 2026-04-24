# On disk representation of BioFormats images in SFE object

\`r lifecycle::badge("experimental")\` At present, the `BioFormatsImage`
is designed for OME-TIFF from Xenium and has not been tested on other
formats that can be read with `BioFormats`. The image is not loaded into
memory, and when it is, the the `BioFormatsImage` object is converted
into
[`ExtImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ExtImage.md)
because the loaded image is of a class that inherits from
[`Image`](https://rdrr.io/pkg/EBImage/man/Image.html). The
[`ExtImage`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ExtImage.md)
class is a thin wrapper inheriting from `VirtualSpatialImage` so it's
compatible with `SpatialExperiment` from which SFE is derived. This
class might drastically change as it matures, say to accommodate other
formats supported by `BioFormats` and to store the transformation matrix
rather than loading image into memory upon transform.

## Usage

``` r
# S4 method for class 'BioFormatsImage'
show(object)

BioFormatsImage(
  path,
  ext = NULL,
  isFull = TRUE,
  origin = c(0, 0),
  transformation = list()
)
```

## Arguments

- object:

  A `BioFormatsImage` object.

- path:

  Path to an OME-TIFF image file.

- ext:

  Numeric vector with names "xmin", "xmax", "ymin", "ymax" in microns
  indicating the spatial extent covered by the image. If `NULL`, then
  the extent will be inferred from the metadata, from physical pixel
  size and the number of pixels.

- isFull:

  Logical, if the extent specified in `ext` is the full extent. If
  `ext = NULL` so it will be inferred from metadata then `isFull = TRUE`
  will be set internally.

- origin:

  Origin of the whole image in the x-y plane, defaults to `c(0,0)`. This
  is shifted when the image is translated. This is not the same as xmin
  and xmax. For example, when the extent is only part of the whole image
  and the whole image itself can be spatially translated, the origin is
  needed to determine which part of the whole image this extent
  corresponds to.

- transformation:

  Named list specifying affine transformation. The list can have names
  "name" and named parameter of the transformation, e.g.
  `list(name = "mirror", direction = "vertical")`, "rotate" and degrees
  = 90 (clockwise), and "scale" and factor = 2. The list can also have
  names "M" for a 2x2 linear transformation matrix in the xy plane and
  "v" for a translation vector of length 2 to specify general affine
  transformation.

## Value

A `BioFormatsImage` object.

## Details

Spatial extent is inferred from OME-TIFF metadata if not specified.
Physical pixel size from the metadata is used to make the extent in
micron space. If physical pixel size is absent from metadata, then the
extent will be in pixel space, which might mean that the image will not
align with the geometries because often the geometry coordinates are in
microns, so a warning is issued in this case.

Affine transformations can be specified in the `transformation`
argument, either by name or by directly specifying the matrix. The
transformations specified by name will always preserve the center of the
image. When named transformations are chained, name and parameter will
be converted to matrix and translation vector the second time a
transformation is specified. If the subsequent transformation happens to
restore the image to its original place, then transformation
specifications will be removed.

## See also

\[isFull()\], \[origin()\]
