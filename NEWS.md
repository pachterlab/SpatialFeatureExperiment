# Version 1.1.6 (04/20/2023)
* Read images as SpatRaster, in read10xVisiumSFE
* read10xVisiumSFE can also convert full resolution image pixels to microns based
on Visium spot spacing
* read10xVisiumSFE no longer transposes output from read10xVisium so the spots
would match the image by default, and to be consistent with SpatialExperiment
* Read standard Vizgen MERFISH output with readVizgen
* SpatRasterImage class inheriting from VirtualSpatialImage for SpatialExperiment compatibility
* Methods of addImg, mirrorImg, and transposeImg for SpatRasterImage and SFE
* Mirror and transpose SFE objects, operating on both geometries and images
* Images are cropped when the SFE object is cropped
* Images are also shifted when removeEmptySpace is called

# Version 1.1.4 (03/02/2023)
* Store SFE package version in object and added SFE method of updateObject to pave way for a potential reimplementation of spatialGraphs.

# Version 1.1.3 (12/20/2022)
* Use BiocNeighbors for k nearest neighbors and distance based neighbors, preserving distance info to avoid slow step to refind distances with sf as done in spdep.
* Added swap_rownames argument in localResult(s) getters so gene symbols from any rowData column can be used to get local results stored under Ensembl IDs.

# Version 1.0.3 (01/11/2023)
* Correctly move the geometries when there are multiple samples
* Use translate = FALSE when using localResult setter for geometries
* More helpful error messages when geometries, localResult, or spatial graphs are absent

# Version 1.0.2 (12/03/2022)
* Correctly move spatialCoords in removeEmptySpace
* Preserve rownames when setting colGeometry for some of all samples

# Version 1.0.0 (11/02/2022)
* First version on Bioconductor

# Version 0.99.4 (09/07/2022)

* Added \code{localResults} field
* Also reimplemented some of the internals behind \code{dimGeometries}

# Version 0.99.0 (02/09/2022)

* Hello world!
* For my personal record, this package was submitted to Bioconductor on July 22
