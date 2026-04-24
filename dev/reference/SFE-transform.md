# Affine transfortaion of SFE object in histological space

These functions perform affine transformations on SFE objects, including
all geometries and images. The transformation is performed on each
sample as a whole. This differs from functions such as
[`mirrorImg`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/mirrorImg.md)
in that `mirrorImg` and `rotateImg` transform the image with the center
at the center of the image itself. In contrast, the center of
transformation here is the center of the bounding box of the entire
sample, including images.

## Usage

``` r
transpose(sfe, sample_id = "all", maxcell = NULL, filename = "")

mirror(
  sfe,
  sample_id = "all",
  direction = c("vertical", "horizontal"),
  maxcell = NULL,
  filename = ""
)

rotate(sfe, sample_id = "all", degrees, maxcell = 1e+07)

translate(sfe, sample_id = "all", v)

scale(sfe, sample_id = "all", factor)

affine(sfe, sample_id = "all", M, v, maxcell = 1e+07)
```

## Arguments

- sfe:

  An SFE object.

- sample_id:

  Sample(s) to transform.

- maxcell:

  Rotating `SpatRasterImage` will convert it into `ExtImage`, loading
  the image into memory. This argument specifies the maximum number of
  pixels in the image loaded into memory. The image will be down sampled
  to approximately this number of pixels.

- filename:

  character. Output filename

- direction:

  character. Should (partially) match "vertical" to flip by rows, or
  "horizontal" to flip by columns

- degrees:

  How many degrees to rotate. Positive number means clockwise and
  negative number means counterclockwise.

- v:

  Vector to spatially translate the SFE object.

- factor:

  Numeric, scaling factor.

- M:

  A 2x2 numeric matrix for the linear transformation in the xy plane.

## Value

An SFE object with the sample(s) transformed.

## Details

For images that are not loaded into memory, `rotateImg` will load
`SpatRasterImage` into memory and all image operations except translate
will load `BioFormatsImage` into memory.

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache
sfe2 <- SpatialFeatureExperiment::transpose(sfe)
sfe3 <- SpatialFeatureExperiment::mirror(sfe)
```
