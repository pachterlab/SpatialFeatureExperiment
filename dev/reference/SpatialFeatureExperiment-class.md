# The SpatialFeatureExperiment class

This class inherits from the
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
(SPE) class, which in turn inherits from
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
(SCE). `SpatialFeatureExperiment` stores geometries of spots or cells in
`sf` objects which form columns of a `DataFrame` which is in turn a
column of the `int_colData` `DataFrame` of the underlying SCE object,
just like `reducedDim` in SCE. Geometries of the tissue outline,
pathologist annotations, and objects (e.g. nuclei segmentation in a
Visium dataset) are stored in `sf` objects in a named list called
`annotGeometries` in `int_metadata`.
