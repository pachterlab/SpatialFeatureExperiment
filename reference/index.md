# Package index

## SpatialFeatureExperiment class

- [`SpatialFeatureExperiment-class`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-class.md)
  : The SpatialFeatureExperiment class
- [`toSpatialFeatureExperiment(`*`<SpatialExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-coercion.md)
  [`toSpatialFeatureExperiment(`*`<SingleCellExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-coercion.md)
  [`toSpatialFeatureExperiment(`*`<Seurat>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-coercion.md)
  : SpatialFeatureExperiment coercion methods
- [`SpatialFeatureExperiment()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment.md)
  : Constructor of SpatialFeatureExperiment object
- [`show(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/show-SpatialFeatureExperiment-method.md)
  : Print method for SpatialFeatureExperiment
- [`unit(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/unit-SpatialFeatureExperiment-method.md)
  : Get unit of a SpatialFeatureExperiment
- [`updateObject(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/updateObject.md)
  [`SFEVersion()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/updateObject.md)
  : Update a SpatialFeatureExperiment object

## Read data into SFE

- [`read10xVisiumSFE()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/read10xVisiumSFE.md)
  : Read 10X Visium data as SpatialFeatureExperiment
- [`readCosMX()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readCosMX.md)
  : Read CosMX data into SFE
- [`readVisiumHD()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readVisiumHD.md)
  : Read Visium HD data
- [`readVizgen()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readVizgen.md)
  : Read Vizgen MERFISH output as SpatialFeatureExperiment
- [`readXenium()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readXenium.md)
  : Read 10X Xenium output as SpatialFeatureExperiment

## Getters and setters

- [`annotGeometries(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`` `annotGeometries<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`annotGeometryNames(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`` `annotGeometryNames<-`( ``*`<SpatialFeatureExperiment>`*`,`*`<character>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`annotGeometry(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`` `annotGeometry<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`tissueBoundary()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  [`` `tissueBoundary<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotGeometries.md)
  : Annotation geometry methods
- [`colFeatureData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colFeatureData.md)
  [`rowFeatureData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colFeatureData.md)
  [`geometryFeatureData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colFeatureData.md)
  [`reducedDimFeatureData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colFeatureData.md)
  : Get global spatial analysis results and metadata of colData,
  rowData, and geometries
- [`colGeometry()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `colGeometry<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`colGeometries()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `colGeometries<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`colGeometryNames()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `colGeometryNames<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`spotPoly()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `spotPoly<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`centroids()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `centroids<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`ROIPoly()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `ROIPoly<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`cellSeg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `cellSeg<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`nucSeg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  [`` `nucSeg<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/colGeometries.md)
  : Column geometry getters and setters
- [`dimGeometries(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dimGeometries.md)
  [`` `dimGeometries<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dimGeometries.md)
  [`dimGeometryNames(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dimGeometries.md)
  [`` `dimGeometryNames<-`( ``*`<SpatialFeatureExperiment>`*`,`*`<numeric>`*`,`*`<character>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dimGeometries.md)
  [`dimGeometry(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dimGeometries.md)
  [`` `dimGeometry<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dimGeometries.md)
  : Dimension geometry methods
- [`getParams()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getParams.md)
  : Get parameters used in spatial methods
- [`localResults(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`` `localResults<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`localResultNames(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`` `localResultNames<-`( ``*`<SpatialFeatureExperiment>`*`,`*`<character>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`localResultFeatures(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`localResultAttrs(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`localResult(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  [`` `localResult<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/localResults.md)
  : Get and set results from local spatial statistics
- [`rowGeometry()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`` `rowGeometry<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`rowGeometries()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`` `rowGeometries<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`rowGeometryNames()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`` `rowGeometryNames<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`txSpots()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  [`` `txSpots<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rowGeometries.md)
  : Row geometry getters and setters
- [`spatialGraphs(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`colGraphs()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`rowGraphs()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`annotGraphs()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `spatialGraphs<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `colGraphs<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `rowGraphs<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `annotGraphs<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`spatialGraphNames(`*`<SpatialFeatureExperiment>`*`,`*`<numeric>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `spatialGraphNames<-`( ``*`<SpatialFeatureExperiment>`*`,`*`<numeric>`*`,`*`<ANY>`*`,`*`<character>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`colGraphNames()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`rowGraphNames()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`annotGraphNames()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `colGraphNames<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `rowGraphNames<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `annotGraphNames<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`spatialGraph(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`colGraph()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`rowGraph()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`annotGraph()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `spatialGraph<-`( ``*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `colGraph<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `rowGraph<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  [`` `annotGraph<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/spatialGraphs.md)
  : Spatial graph methods

## Preprocessing and QC

- [`findDebrisCells(`*`<matrix>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findDebrisCells.md)
  [`findDebrisCells(`*`<sfc>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findDebrisCells.md)
  [`findDebrisCells(`*`<sf>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findDebrisCells.md)
  [`findDebrisCells(`*`<SpatialExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findDebrisCells.md)
  : Identify cells in small bits outside the main piece of tissue
- [`getTissueBoundaryConcave(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryConcave.md)
  [`getTissueBoundaryConcave(`*`<sfc>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryConcave.md)
  [`getTissueBoundaryConcave(`*`<sf>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryConcave.md)
  : Get tissue boundary from concave hull of cell geometries
- [`getTissueBoundaryImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryImg.md)
  : Get tissue boundary from histology image

## Find spatial neighborhood graph

- [`findSpatialNeighbors(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findSpatialNeighbors.md)
  : Find spatial neighborhood graph
- [`findVisiumGraph()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findVisiumGraph.md)
  : Find spatial neighborhood graphs for Visium spots
- [`findVisiumHDGraph()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/findVisiumHDGraph.md)
  : Find Visium HD spatial neighborhood graph

## Transcript spots

- [`formatTxSpots()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/formatTxSpots.md)
  [`addTxSpots()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/formatTxSpots.md)
  : Read and process transcript spots geometry for SFE
- [`formatTxTech()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/formatTxTech.md)
  [`addTxTech()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/formatTxTech.md)
  : Read and process transcript spots for specific commercial
  technologies
- [`readSelectTx()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readSelectTx.md)
  [`addSelectTx()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/readSelectTx.md)
  : Read transcript spots of select genes

## Geometric/spatial operations

- [`transpose()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-transform.md)
  [`mirror()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-transform.md)
  [`rotate()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-transform.md)
  [`translate()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-transform.md)
  [`scale()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-transform.md)
  [`affine()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-transform.md)
  : Affine transfortaion of SFE object in histological space
- [`addVisiumSpotPoly()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/addVisiumSpotPoly.md)
  : Add Visium spot polygons to colGeometry
- [`aggregate(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/aggregate-SpatialFeatureExperiment-method.md)
  : Aggregate data in SFE using geometry
- [`aggregateTx()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/aggregateTx.md)
  [`aggregateTxTech()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/aggregateTx.md)
  : Aggregate transcript spots from file
- [`annotOp()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotOp.md)
  : Binary operations for geometry of each cell/spot and annotation
- [`annotPred()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotPred.md)
  [`annotNPred()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotPred.md)
  : Binary predicates for geometry of each cell/spot and annotation
- [`annotSummary()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/annotSummary.md)
  : Summarize attributes of an annotGeometry for each cell/spot
- [`bbox(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/bbox-SpatialFeatureExperiment-method.md)
  : Find bounding box of SFE objects
- [`crop()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/crop.md)
  : Crop an SFE object with a geometry
- [`removeEmptySpace()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/removeEmptySpace.md)
  : Remove empty space
- [`rotateMinRect(`*`<matrix>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateMinRect.md)
  [`rotateMinRect(`*`<sfc>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateMinRect.md)
  [`rotateMinRect(`*`<sf>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateMinRect.md)
  [`rotateMinRect(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateMinRect.md)
  : Rotate the tissue to reduce empty space in plots
- [`splitByCol(`*`<SpatialFeatureExperiment>`*`,`*`<sf>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitByCol.md)
  [`splitByCol(`*`<SpatialFeatureExperiment>`*`,`*`<sfc>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitByCol.md)
  [`splitByCol(`*`<SpatialFeatureExperiment>`*`,`*`<list>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitByCol.md)
  [`splitSamples()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitByCol.md)
  [`splitContiguity()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitByCol.md)
  : Split SFE object with categorical vector or geometry
- [`splitComponent(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitComponent.md)
  [`splitComponent(`*`<sf>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitComponent.md)
  [`splitComponent(`*`<sfc>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/splitComponent.md)
  : Split by graph component
- [`st_any_pred()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/st_any_pred.md)
  [`st_any_intersects()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/st_any_pred.md)
  [`st_n_pred()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/st_any_pred.md)
  [`st_n_intersects()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/st_any_pred.md)
  : Simple geometry predicates

## Non-spatial operations

- [`` `[`( ``*`<SpatialFeatureExperiment>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatialFeatureExperiment-subset.md)
  : Subsetting SpatialFeatureExperiment objects
- [`cbind(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cbind-SpatialFeatureExperiment-method.md)
  : Concatenate SpatialFeatureExperiment objects

## Image classes

- [`isFull(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/BioFormatsImage-getters.md)
  [`origin(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/BioFormatsImage-getters.md)
  [`transformation(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/BioFormatsImage-getters.md)
  :

  Other `BioFormatsImage` getters

- [`show(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/BioFormatsImage.md)
  [`BioFormatsImage()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/BioFormatsImage.md)
  : On disk representation of BioFormats images in SFE object

- [`show(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ExtImage.md)
  [`ExtImage()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ExtImage.md)
  :

  Use the EBImage `Image` class in SFE objects

- [`SpatRasterImage()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatRasterImage.md)
  [`show(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SpatRasterImage.md)
  : SpatRaster representation of images in SFE objects

- [`toExtImage(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/toExtImage.md)
  [`toExtImage(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/toExtImage.md)
  : Convert images to ExtImage

- [`toSpatRasterImage(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/toSpatRasterImage.md)
  [`toSpatRasterImage(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/toSpatRasterImage.md)
  : Convert images to SpatRasterImage

## Image methods

- [`` `Img<-`( ``*`<SpatialExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/Img-set-SpatialExperiment-method.md)
  : Image setter
- [`addImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  [`transposeImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  [`mirrorImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  [`rotateImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  [`translateImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  [`scaleImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  [`affineImg(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/SFE-image.md)
  : Methods for handling image-related data
- [`cropImg(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md)
  [`cropImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md)
  [`cropImg(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/cropImg.md)
  : Crop images
- [`dim(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-BioFormatsImage-method.md)
  : Find dimension of BioFormatsImage
- [`dim(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/dim-ExtImage-method.md)
  : Find dimensions of ExtImage
- [`ext(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
  [`ext(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
  [`ext(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
  [`` `ext<-`( ``*`<BioFormatsImage>`*`,`*`<numeric>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
  [`` `ext<-`( ``*`<ExtImage>`*`,`*`<numeric>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
  [`` `ext<-`( ``*`<SpatRasterImage>`*`,`*`<numeric>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/ext.md)
  : Get and set extent of image objects
- [`imgRaster(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md)
  [`imgRaster(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md)
  [`imgRaster(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgRaster.md)
  : Get the image from \*Image class
- [`imgSource(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md)
  [`imgSource(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md)
  [`imgSource(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imgSource.md)
  : Source of images that are on disk

## Image affine transformations

- [`affineImg(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md)
  [`affineImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md)
  [`affineImg(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/affineImg.md)
  : Affine transformation of images
- [`mirrorImg(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md)
  [`mirrorImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md)
  [`mirrorImg(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/mirrorImg.md)
  : Mirror/flip images
- [`rotateImg(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md)
  [`rotateImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md)
  [`rotateImg(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/rotateImg.md)
  : Rotate image
- [`scaleImg(`*`<AlignedSpatialImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md)
  [`scaleImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/scaleImg.md)
  : Scale image
- [`translateImg(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md)
  [`translateImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md)
  [`translateImg(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/translateImg.md)
  : Translate/shift image in space
- [`transposeImg(`*`<SpatRasterImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
  [`transposeImg(`*`<BioFormatsImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
  [`transposeImg(`*`<ExtImage>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/transposeImg.md)
  : Transpose images

## Utilities

- [`aggBboxes()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/aggBboxes.md)
  : Aggregate bounding boxes
- [`bbox_center()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/bbox_center.md)
  : Find center of bounding box
- [`changeSampleIDs()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/changeSampleIDs.md)
  : Change sample IDs
- [`containsOutOfMemoryData(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/containsOutOfMemoryData-SpatialFeatureExperiment-method.md)
  : Whether an SFE object contains out of memory data
- [`df2sf()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md)
  : From ordinary data frame to sf to construct SFE object
- [`gdalParquetAvailable()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/gdalParquetAvailable.md)
  : Check if Parquet GDAL driver is available
- [`getPixelSize()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getPixelSize.md)
  : Get physical size of pixels
- [`getTechTxFields()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTechTxFields.md)
  : Get relevant fields and file paths for transcript spots
- [`imageIDs()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/imageIDs.md)
  : Show all image_ids in the SFE object
- [`multi_listw2sparse()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/multi_listw2sparse.md)
  : Convert multiple listw graphs into a single sparse adjacency matrix
- [`sampleIDs()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/sampleIDs.md)
  : Get all unique sample IDs
- [`saveRDS(`*`<SpatialFeatureExperiment>`*`)`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/saveRDS-SpatialFeatureExperiment-method.md)
  : Save SpatialFeatureExperiment as RDS file
- [`visium_row_col`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/visium_row_col.md)
  : Row and columns of Visium barcodes on the slide

## Re-exports

- [`colData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`rowData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`` `colData<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`spatialCoords()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`` `spatialCoords<-`() ``](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`spatialCoordsNames()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`getImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`imgData()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`rmvImg()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`counts()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`logcounts()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  [`reducedDim()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/reexports.md)
  : Functions re-exported from other packages
