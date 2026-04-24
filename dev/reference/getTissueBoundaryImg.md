# Get tissue boundary from histology image

This function gets the tissue boundary from image and makes sure that it
is properly aligned with the geometries in the SFE object. Note that it
will not work on fluorescent images when the fluorescence looks sparse.
In that case, use
[`getTissueBoundaryConcave`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/getTissueBoundaryConcave.md).

## Usage

``` r
getTissueBoundaryImg(
  sfe,
  sample_id = NULL,
  image_id = NULL,
  image_type = c("brightfield", "fluorescent"),
  channel = NULL,
  n_pieces = 1,
  resolution = 4,
  maxcell = 1e+07,
  fill_holes = FALSE,
  simplify = TRUE,
  dTolerance = 0
)
```

## Arguments

- sfe:

  An SFE object with images

- sample_id:

  Sample id(s) whose tissue boundaries are to be found.

- image_id:

  ID of image to use to get boundary.

- image_type:

  Character, either "brightfield" or "fluorescent"

- channel:

  Channel to use for tissue segmentation. If `NULL` use average of all
  channels.

- n_pieces:

  Number of pieces; only this number of largest pieces are kept. Smaller
  pieces will be considered debris and removed. Can be a vector of the
  same length as `sample_id`, or if `sample_id` is "all" then same
  length as the number of samples in the SFE object to specify different
  number of pieces in different samples. If `n_pieces` is length 1 while
  there are multiple samples, then the same number is applied to all
  samples.

- resolution:

  Integer, which resolution to use for tissue boundary in a pyramidal
  OME-TIFF stack. Only applicable to `BioFormatsImage`. Note that the
  image will be loaded into memory and you usually don't need the
  highest resolution for the tissue boundary.

- maxcell:

  Max number of pixels when loading `SpatRasterImage` into memory.

- fill_holes:

  Logical, whether to fill holes in the tissue, to only get the outer
  outline.

- simplify:

  Logical, whether to simplify the output polygon.

- dTolerance:

  Distance tolerance when simplifying the polygon, in the same unit as
  the geometries in the SFE object.

## Value

A `sf` data frame with columns `sample_id` and `geometry`.
