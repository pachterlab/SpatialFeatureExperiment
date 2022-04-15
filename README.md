
# SpatialFeatureExperiment

<!-- badges: start -->
<!-- badges: end -->

`SpatialFeatureExperiment` (SFE) is a new S4 class built on top of `SpatialExperiment` (SPE) that incorporates geometries and geometry operations with the `sf` package. Examples of such geometries are Visium spots represented with polygons representing their actual size, cell or nuclei segmentation polygons, tissue boundary polygons, pathologist annotation of histological regions, and transcript spots of genes. Using `sf`, `SpatialFeatureExpzeriment` leverages well-established and optimized C++ libraries underlying `sf` for geometry operations such as finding whether geometries intersect, finding the intersection geometries, buffering a geometry with a margin, and etc. A schematic of the SFE object is shown below:

<img src="vignettes/sfe_schematics.png" width="800"/>

These are the additions to the SPE object:

* `colGeometries` are `sf` data frames associated with the entities that correspond to columns of the gene count matrix, such as Visium spots or cells. The geometries in the `sf` data frames can be Visium spot centroids, Visium spot polygons, or for dataset with single cell resolution, cell or nuclei segmentations. Multiple `colGeometries` can be stored in the same SFE object, such as one for cell segmentation and another for nuclei segmentation. There can be non-spatial, attribute columns in a `colGeometry` rather than `colData`, because the `sf` class allows users to specify how attributes relate to geometries, such as "constant", "aggregate", and "identity". See the `agr` argument of the [`st_sf` documentation](https://r-spatial.github.io/sf/reference/sf.html).
* `colGraphs` are spatial neighborhood graphs of cells or spots. The graphs have class `listw` (`spdep` package) and the `colPairs` of `SingleCellExperiment` was not used so no conversion is necessary to use the numerous spatial dependency functions from `spdep`, such as those for Moran's I, Geary's C, Getis-Ord Gi*, LOSH, and etc., as well as other classical spatial statistics packages such as `spatialreg` and `adespatial`.
* `rowGeometries` are similar to `colGeometries`, but are for entities that correspond to rows of the gene count matrix, such as genes. A potential use case would be transcript spots for each gene in smFISH or in situ sequencing based datasets.
* `rowGraphs` are similar to `colGraphs`. A potential use case may be spatial colocalization of transcripts of different genes.
* `annotGeometries` are `sf` data frames associated with the dataset but not directly with the gene count matrix, such as tissue boundaries, histological regions, cell or nuclei segmentation in Visium datasets, and etc. These geometries are stored in this object to facilitate plotting and using `sf` for operations such as to find the number of nuclei in each Visium spot and which histological regions each Visium spot intersects. Unlike `colGeometries` and `rowGeometries`, the number of rows in the `sf` data frames in `annotGeometries` is not constrained by the dimension of the gene count matrix and can be arbitrary.
* `annotGraphs` are similar to `colGraphs` and `rowGraphs`, but are for entities not directly associated with the gene count matrix, such as spatial neighborhood graphs for nuclei in a Visium dataset or other objects like myofibers. These graphs are relevant to `spdep` analyses of attributes of these geometries such as spatial autocorrelation in morphological metrics of myofibers and nuclei. With geometry operations with `sf`, these attributes and results of analyses of these attributes (e.g. spatial regions defined by the attributes) may be related back to gene expression.

## Installation

You can install the development version of SpatialFeatureExperiment from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pachterlab/SpatialFeatureExperiment")
```

This package will be submitted to Bioconductor after we write examples and vignettes and polish the functions after using it to perform more analyses.
