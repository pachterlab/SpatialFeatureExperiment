
# SpatialFeatureExperiment

<!-- badges: start -->
<!-- badges: end -->

`SpatialFeatureExperiment` is a new S4 class built on top of `SpatialExperiment` that incorporates geometries and geometry operations with the `sf` package. Examples of such geometries are Visium spots represented with polygons representing their actual size, cell or nuclei segmentation polygons, tissue boundary polygons, pathologist annotation of histological regions, and transcript spots of genes. Using `sf`, `SpatialFeatureExpzeriment` leverages well-established and optimized C++ libraries underlying `sf` for geometry operations such as finding whether geometries intersect, finding the intersection geometries, buffering a geometry with a margin, and etc.

## Installation

You can install the development version of SpatialFeatureExperiment from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pachterlab/SpatialFeatureExperiment")
```

This package will be submitted to Bioconductor after I write unit tests, examples, and vignettes.
