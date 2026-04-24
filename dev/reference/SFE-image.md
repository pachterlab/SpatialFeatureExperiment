# Methods for handling image-related data

Generics of these functions are defined in `SpatialExperiment`, except
for `transposeImg`. These SFE methods cater to the new image-related
classes in SFE. The SPE method for `getImg`, `rmvImg`, and `imgRaster`
don't need to be modified for SFE and are hence not implemented here,
but are simply re-exported.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
addImg(x, imageSource, sample_id = 1L, image_id, extent = NULL, scale_fct = 1)

# S4 method for class 'SpatialFeatureExperiment'
transposeImg(
  x,
  sample_id = 1L,
  image_id = NULL,
  maxcell = 1e+07,
  filename = ""
)

# S4 method for class 'SpatialFeatureExperiment'
mirrorImg(
  x,
  sample_id = 1L,
  image_id = NULL,
  direction = "vertical",
  maxcell = 1e+07,
  filename = ""
)

# S4 method for class 'SpatialFeatureExperiment'
rotateImg(x, sample_id = 1L, image_id = NULL, degrees, maxcell = 1e+07)

# S4 method for class 'SpatialFeatureExperiment'
translateImg(x, sample_id = 1L, image_id = NULL, v)

# S4 method for class 'SpatialFeatureExperiment'
scaleImg(x, sample_id = 1L, image_id = NULL, factor)

# S4 method for class 'SpatialFeatureExperiment'
affineImg(x, sample_id = 1L, image_id = NULL, M, v)
```

## Arguments

- x:

  A SFE object.

- imageSource:

  a character string specifying an image file name (.png, .jpg or .tif)
  or URL to source the image from

- sample_id:

  Which sample the image is associated with. Use
  [`sampleIDs`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/sampleIDs.md)
  to get sample IDs present in the SFE object.

- image_id:

  Image ID, such as "lowres" and "hires" for Visium data and "DAPI" and
  "PolyT" for Vizgen MERFISH data.

- extent:

  A numeric vector of length 4 with names of the set xmin, ymin, xmax,
  and ymax, specifying the extent of the image.

- scale_fct:

  Scale factor – multiply pixel coordinates in full resolution image by
  this scale factor should yield pixel coordinates in a different
  resolution. `extent` takes precedence over `scale_fct`.

- maxcell:

  Max number of pixels to load `SpatRasterImage` into memory. The
  default 1e7 is chosen because this is the approximate number of pixels
  in the medium resolution image at `resolution = 4L` in Xenium OME-TIFF
  to make different methods of this function consistent.

- filename:

  character. Output filename

- direction:

  character. Should (partially) match "vertical" to flip by rows, or
  "horizontal" to flip by columns

- degrees:

  How many degrees to rotate. Positive number means clockwise and
  negative number means counterclockwise.

- v:

  A numeric vector of length 2 specifying the vector in the xy plane to
  translate the SFE object.

- factor:

  Numeric, scaling factor.

- M:

  A 2x2 numeric matrix for the linear transformation in the xy plane.

## Details

Method of
[`transposeImg`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/transposeImg.md),
[`mirrorImg`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/mirrorImg.md),
and
[`rotateImg`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/rotateImg.md)
perform the method on all images within the SFE object that are
specified with `sample_id` and `image_id`. For images that are not
loaded into memory, `rotateImg` will load `SpatRasterImage` into memory
and all image operations except translate will load `BioFormatsImage`
into memory.

## Note

If the image is already a GeoTIFF file that already has an extent, then
the extent associated with the file will be honored and the `extent` and
`scale_fct` arguments are ignored. Transposing the image is just like
transposing a matrix. It's flipped about the line going from the top
left to the bottom right.

## See also

Other image methods:
[`affineImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/affineImg.md),
[`cropImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/cropImg.md),
[`dim,BioFormatsImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-BioFormatsImage-method.md),
[`dim,ExtImage-method`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/dim-ExtImage-method.md),
[`ext()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/ext.md),
[`imgRaster()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/imgRaster.md),
[`imgSource()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/imgSource.md),
[`mirrorImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/mirrorImg.md),
[`rotateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/rotateImg.md),
[`scaleImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/scaleImg.md),
[`translateImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/translateImg.md),
[`transposeImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/dev/reference/transposeImg.md)

## Examples

``` r
library(SFEData)
sfe <- McKellarMuscleData("small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
img_path <- system.file(file.path("extdata", "sample01", "outs", "spatial",
"tissue_lowres_image.png"), package = "SpatialFeatureExperiment")
sfe <- addImg(sfe, img_path, sample_id = "Vis5A", image_id = "lowres", scale_fct =
0.023)
img <- getImg(sfe)
# SpatRasterImage method
img_t <- transposeImg(img)
# SFE method
sfe <- transposeImg(sfe, sample_id = "Vis5A", image_id = "lowres")
```
