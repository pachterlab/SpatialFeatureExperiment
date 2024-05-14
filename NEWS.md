# Version 1.6.1 (05/09/2024)
* readRDS converts old style SpatRasterImage to the new style
* readSelectTx and addSelectTx functions to read transcript spots from a few select genes from the parquet output of formatTxSpots or add them to an SFE object
* Added formatTxTech and addTxTech functions, basically thin wrappers of formatTxSpots and addTxSpots with presets for Vizgen, Xenium, and CosMX

# Version 1.6.0 (04/29/2024)
* Changed defaults from sample_id = NULL to sample_id = 1L when dealing with 1 sample or "all" when dealing with multiple samples
* dim method for BioFormatsImage that doesn't load the image into memory
* Deal with univariate spatial results in featureData in cbind and changeSampleID
* Fixed super embarrassing bug in cbind that fails when combining more than 2 SFE objects
* Updated readXenium for XOA v2
* Updated BioFormatsImage to store affine transform info rather than converting to EBImage after transform
* Speed up affine transformation of sf geometries with sfheaders
* Coercion from Seurat to SFE
* SpatRasterImage and EBImage directly inherit from SpatRaster and Image respectively so the user no longer needs to call imgRaster every time they plots or operates on the image, which I find really annoying.
* Changed name EBImage to ExtImage to reduce confusion
* Bug fixes on image affine transformation
* Exporting some util functions: aggBboxes, getPixelSize, and imageIDs
* Read select genes as rowGeometry and add to existing rowGeometry without erasing existing genes in the same rowGeometry

# Version 1.5.2 (03/04/2024)
* Added readXenium (for XOA v1)
* Added BioFormatsImage and EBImage classes to deal with Xenium OME-TIFF
* Conversion between SpatRasterImage, BioFormatsImage, and EBImage
* Overhaul of geometry operation functions for images and SFE objects for the new image classes, including bbox, crop, and affine transforms
* Don't throw error when there are no rows or columns left after [ subsetting
* cbind for multiple samples that have rowGeometry
* Rewrote df2df with the much faster sfheaders, deprecating the less efficient BPPARAM argument

# Version 1.5.1 (02/02/2024)
* Added support for rowGeometry and transcript spots
* Reformat transcript spot files from Vizgen and CosMX
* Improved readVizgen for transcript spots
* Added readCosMX

# Version 1.3.1 (09/22/2023)
* Refactored to remove "missing" methods for geometries, graphs, and local results.
* Changed defaults from sample_id = NULL to sample_id = "all" unless only one sample can be specified.

# Version 1.2.3 (08/18/2023)
* Fixed bug when Visium graph is not added when filtered matrix from only one sample is read with read10xVisiumSFE.
* Changed the way pixels are converted to microns in Visium. Old way: use top left corner of Visium spot array to compute spacing between spots, doesn't work for filtered data when there're singleton spots. New way: Use median row/col indices, more robust when there're singletons. When spacing is used for the conversion, spot size is found to vary across datasets.
* Added saveRDS method for SFE objects, so SpatRaster images are wrapped before saving and unwrapped on the fly when they're requested.
* Fixed bug when the wrong bounding boxes are used to crop images when SFE object is subsetted and there're multiple samples.

# Version 1.2.2 (07/21/2023)
* Fixed embarrassing documentation mismatch in localResults

# Version 1.2.1 (04/26/2023)
* Fixed bug in .check_features and .symbol2id where "symbol" column is hard coded

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
