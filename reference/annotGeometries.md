# Annotation geometry methods

"Annotation geometry" refers to Simple Feature (`sf`) geometries NOT
associated with rows (features, genes) or columns (cells or spots) of
the gene count matrix in the `SpatialFeatureExperiment` object. So there
can be any number of rows in the `sf` data frame specifying the
geometry. Examples of such geometries are tissue boundaries, pathologist
annotation of histological regions, and objects not characterized by
columns of the gene count matrix (e.g. nuclei segmentation in a Visium
dataset where the columns are Visium spots). This page documents getters
and setters for the annotation geometries. Internally, annotation
geometries are stored in `int_metadata`.

## Usage

``` r
# S4 method for class 'SpatialFeatureExperiment'
annotGeometries(x)

# S4 method for class 'SpatialFeatureExperiment'
annotGeometries(x, translate = TRUE, ...) <- value

# S4 method for class 'SpatialFeatureExperiment'
annotGeometryNames(x)

# S4 method for class 'SpatialFeatureExperiment,character'
annotGeometryNames(x) <- value

# S4 method for class 'SpatialFeatureExperiment'
annotGeometry(x, type = 1L, sample_id = NULL)

# S4 method for class 'SpatialFeatureExperiment'
annotGeometry(x, type = 1L, sample_id = NULL, translate = TRUE, ...) <- value

tissueBoundary(x, sample_id = 1L)

tissueBoundary(x, sample_id = 1L, translate = TRUE, ...) <- value
```

## Arguments

- x:

  A `SpatialFeatureExperiment` object.

- translate:

  Logical. Only used if
  [`removeEmptySpace`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/removeEmptySpace.md)
  has been run of the SFE object. If that's the case, this argument
  indicates whether the new value to be assigned to the geometry is in
  the coordinates prior to removal of empty space so it should be
  translated to match the new coordinates after removing empty space.
  Default to `TRUE`.

- ...:

  `spatialCoordsNames, spotDiameter, geometryType` passed to
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md).
  Defaults are the same as in
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md).
  For `dimGeometries<-` only: `geometryType` can be a character vector
  of the geometry type of each data frame in the list of the same length
  as the list if the data frames specify different types of geometries.

- value:

  Value to set. For `annotGeometry`, must be a `sf` data frame, or an
  ordinary data frame that can be converted to a `sf` data frame (see
  [`df2sf`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/df2sf.md)).
  For `annotGeometries`, must be a list of such `sf` or ordinary data
  frames. There must be a column `sample_id` to indicate the sample the
  geometries are for, and the `sample_id` must also appear in `colData`.

- type:

  An integer specifying the index or string specifying the name of the
  \*Geometry to query or replace. If missing, then the first item in the
  \*Geometries will be returned or replaced.

- sample_id:

  Sample ID to get or set geometries.

## Value

Getters for multiple geometries return a named list. Getters for names
return a character vector of the names. Getters for single geometries
return an `sf` data frame. Setters return an SFE object.

## Details

Wrapper for getter and setter of special geometry:

- tisseuBoundary:

  Boundary of the tissue of interest, including holes. This is usually
  of geometry type MULTIPOLYGON, though geometries in `annotGeometries`
  can have any type supported by `sf`.

## Examples

``` r
# Example dataset
library(SFEData)
sfe_small <- McKellarMuscleData(dataset = "small")
#> see ?SFEData and browseVignettes('SFEData') for documentation
#> loading from cache

# Get all annotation geometries, returning a named list
annotGeometries(sfe_small)
#> $tissueBoundary
#> Simple feature collection with 1 feature and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 5094 ymin: 13000 xmax: 7000 ymax: 14969
#> CRS:           NA
#>   ID                       geometry sample_id
#> 7  7 POLYGON ((5094 13000, 5095 ...     Vis5A
#> 
#> $myofiber_full
#> Simple feature collection with 195 features and 2 fields
#> Geometry type: GEOMETRY
#> Dimension:     XY
#> Bounding box:  xmin: 5072 ymin: 13000 xmax: 7000 ymax: 14956
#> CRS:           NA
#> First 10 features:
#>      lyr.1                       geometry sample_id
#> 1018  1018 POLYGON ((6620 14908, 6628 ...     Vis5A
#> 1021  1021 POLYGON ((6476 13416, 6492 ...     Vis5A
#> 1024  1024 GEOMETRYCOLLECTION (LINESTR...     Vis5A
#> 1041  1041 POLYGON ((6024 14592, 6044 ...     Vis5A
#> 1052  1052 POLYGON ((5948 14664, 5976 ...     Vis5A
#> 1058  1058 POLYGON ((6184 13916, 6192 ...     Vis5A
#> 1131  1131 POLYGON ((6560 14040, 6572 ...     Vis5A
#> 1149  1149 POLYGON ((6948 14792, 6956 ...     Vis5A
#> 1164  1164 GEOMETRYCOLLECTION (LINESTR...     Vis5A
#> 1198  1198 POLYGON ((7000 14088, 6988 ...     Vis5A
#> 
#> $myofiber_simplified
#> Simple feature collection with 195 features and 34 fields
#> Geometry type: GEOMETRY
#> Dimension:     XY
#> Bounding box:  xmin: 5075.2 ymin: 13000 xmax: 7000 ymax: 14956
#> CRS:           NA
#> First 10 features:
#>      lyr.1   area perimeter eccentricity       theta  sine_theta convexity
#> 1018  1018 1254.5  144.3111    0.7556634 -1.43691705 -0.99105155 0.9439428
#> 1021  1021 1294.0  154.6381    0.8373427  1.22884082  0.94210073 0.8970537
#> 1024  1024 1071.0  137.5860    0.7616274 -0.56418348 -0.53472602 0.8876917
#> 1041  1041 1117.5  143.9801    0.4602999  0.42465952  0.41201057 0.8830502
#> 1052  1052  702.5  108.8714    0.5921893  0.05022803  0.05020692 0.9117456
#> 1058  1058  748.0  130.9098    0.7333168  1.48690171  0.99648291 0.7877830
#> 1131  1131  921.0  123.1718    0.4938581  0.08929549  0.08917687 0.9109792
#> 1149  1149 1350.0  155.2466    0.8064659  1.45088782  0.99281958 0.9127789
#> 1164  1164 1361.0  160.5250    0.5199006 -1.43088381 -0.99022820 0.8668790
#> 1198  1198 1359.5  147.4640    0.7536808 -1.55943349 -0.99993544 0.9493715
#>        h.asm.s1 h.con.s1  h.cor.s1  h.var.s1  h.idm.s1 h.sav.s1 h.sva.s1
#> 1018 0.04569305 8.823872 0.6175303 12.535387 0.4677583 55.93346 3043.557
#> 1021 0.06676416 2.971181 0.3866366  3.422040 0.5244005 56.36618 3086.581
#> 1024 0.04072365 5.629854 0.5509958  7.269267 0.4544623 53.67451 2789.833
#> 1041 0.03789957 7.404033 0.6839055 12.711739 0.4583396 54.40403 2876.513
#> 1052 0.05586671 5.919123 0.5754451  7.970974 0.5074887 55.10506 2951.151
#> 1058 0.06039099 4.356942 0.6920377  8.073825 0.5117327 56.48317 3107.141
#> 1131 0.05073276 4.830403 0.6019241  7.067188 0.4921981 56.95746 3152.026
#> 1149 0.05089796 5.356243 0.7815510 13.259709 0.4837150 55.53230 3015.888
#> 1164 0.05300985 4.562002 0.6782312  8.088945 0.4944515 55.49526 2993.993
#> 1198 0.04410018 5.977108 0.5276979  7.327632 0.4641892 55.45479 2981.551
#>       h.sen.s1 h.ent.s1 h.dva.s1  h.den.s1   h.f12.s1  h.f13.s1   h.asm.s2
#> 1018 1.1042774 1.702237 8.823872 0.7768725 0.08287254 0.3699324 0.04843700
#> 1021 0.8696623 1.375543 2.971181 0.6110582 0.03857801 0.2295028 0.07371797
#> 1024 1.0400578 1.655837 5.629854 0.7268132 0.05447687 0.2975921 0.04179965
#> 1041 1.1398830 1.746922 7.404033 0.7499876 0.10322951 0.4161292 0.04108058
#> 1052 0.9831099 1.515003 5.919123 0.6924231 0.06404495 0.3088334 0.05065103
#> 1058 0.9565778 1.475250 4.356942 0.6600786 0.09843866 0.3763597 0.06705178
#> 1131 0.9879345 1.523556 4.830403 0.6717902 0.06896841 0.3211188 0.04611804
#> 1149 1.0143578 1.572877 5.356243 0.6964793 0.09582218 0.3826318 0.05029586
#> 1164 0.9957274 1.548318 4.562002 0.6773656 0.08123645 0.3505345 0.05579443
#> 1198 1.0285324 1.638017 5.977108 0.7321981 0.06439037 0.3213296 0.04349491
#>       h.con.s2  h.cor.s2  h.var.s2  h.idm.s2 h.sav.s2 h.sva.s2  h.sen.s2
#> 1018 13.956859 0.3591553 11.889424 0.4635219 56.08628 3053.482 1.0964515
#> 1021  4.374046 0.1620524  3.609976 0.5356106 56.44869 3096.160 0.8600733
#> 1024  7.633299 0.2524210  6.105346 0.4313602 53.81512 2802.705 0.9954560
#> 1041  9.969246 0.4542513 10.133550 0.4460751 54.74504 2904.280 1.1007054
#> 1052  8.709208 0.3130783  7.339302 0.4888779 55.15509 2949.709 0.9973425
#> 1058  5.844347 0.4605127  6.416574 0.5263770 56.37004 3088.457 0.9386880
#> 1131  6.231405 0.3613120  5.878286 0.4723205 57.09327 3162.454 0.9757856
#> 1149  9.113317 0.6047476 12.528476 0.4714423 55.47891 3003.012 1.0178087
#> 1164  5.682109 0.4433550  6.103889 0.4904749 55.70128 3011.191 0.9613706
#> 1198  8.985191 0.2719318  7.170570 0.4435528 55.39049 2970.667 1.0308582
#>      h.ent.s2  h.dva.s2  h.den.s2   h.f12.s2  h.f13.s2
#> 1018 1.648689 13.956859 0.8096457 0.07724304 0.3522316
#> 1021 1.347868  4.374046 0.6222874 0.02909852 0.1975308
#> 1024 1.592808  7.633299 0.7431862 0.03346647 0.2297161
#> 1041 1.688399  9.969246 0.7857339 0.09662418 0.3969055
#> 1052 1.501903  8.709208 0.6950301 0.05405809 0.2829466
#> 1058 1.423354  5.844347 0.6871059 0.07950027 0.3334109
#> 1131 1.503655  6.231405 0.7083074 0.05990168 0.2977786
#> 1149 1.553251  9.113317 0.7371693 0.08636114 0.3616774
#> 1164 1.496062  5.682109 0.6900179 0.06215469 0.3025037
#> 1198 1.645228  8.985191 0.7904764 0.05511628 0.2983840
#>                            geometry sample_id
#> 1018 POLYGON ((6620 14908, 6628 ...     Vis5A
#> 1021 POLYGON ((6476 13416, 6492 ...     Vis5A
#> 1024 POLYGON ((6980 13492, 6992 ...     Vis5A
#> 1041 POLYGON ((6024 14592, 6044 ...     Vis5A
#> 1052 POLYGON ((6004 14652, 6016 ...     Vis5A
#> 1058 POLYGON ((6184 13916, 6192 ...     Vis5A
#> 1131 POLYGON ((6560 14040, 6592 ...     Vis5A
#> 1149 POLYGON ((6948 14792, 6956 ...     Vis5A
#> 1164 POLYGON ((7000 13152, 6984 ...     Vis5A
#> 1198 POLYGON ((7000 14088, 6988 ...     Vis5A
#> 
#> $nuclei
#> Simple feature collection with 442 features and 9 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 5097.075 ymin: 13000 xmax: 7000 ymax: 14948.44
#> CRS:           NA
#> First 10 features:
#>          ID index     area roundness eccentricity aspect_ratio     angle
#> 10084 10084 10084 238.4863  1.085386    0.7328259     1.469689 139.27116
#> 10175 10175 10175 493.0881  1.067844    0.6938556     1.388669 148.65620
#> 10181 10181 10181 177.5525  1.061164    0.6553132     1.323877 122.36570
#> 10189 10189 10189 323.6471  1.153316    0.8276987     1.781982  45.38620
#> 10246 10246 10246 416.4098  1.072617    0.7292323     1.461420  97.54036
#> 10316 10316 10316 269.0076  1.067543    0.7193451     1.439568  11.75157
#> 10426 10426 10426 375.6754  1.028084    0.1588077     1.012854  16.43842
#> 10432 10432 10432 225.9670  1.038551    0.4909133     1.147831 162.73535
#> 10440 10440 10440 310.9183  1.052642    0.6360212     1.295888  72.76263
#> 10486 10486 10486 329.1423  1.084983    0.7364244     1.478153 159.36537
#>       convexity                       geometry sample_id
#> 10084 0.9920712 POLYGON ((5776.277 14250, 5...     Vis5A
#> 10175 0.9944043 POLYGON ((6920.169 14910, 6...     Vis5A
#> 10181 0.9914434 POLYGON ((6414.457 14760, 6...     Vis5A
#> 10189 0.9869685 POLYGON ((6396.188 13995, 6...     Vis5A
#> 10246 0.9953556 POLYGON ((5838.198 14742, 5...     Vis5A
#> 10316 0.9958743 POLYGON ((5407.794 13000, 5...     Vis5A
#> 10426 0.9945204 POLYGON ((6014.04 14542, 60...     Vis5A
#> 10432 0.9912810 POLYGON ((6308.317 14453, 6...     Vis5A
#> 10440 0.9936670 POLYGON ((6145.372 13068, 6...     Vis5A
#> 10486 0.9916191 POLYGON ((6091.856 13626, 6...     Vis5A
#> 
#> $nuclei_centroid
#> Simple feature collection with 434 features and 8 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 5103 ymin: 13004 xmax: 7000 ymax: 14942
#> CRS:           NA
#> First 10 features:
#>       index     area roundness eccentricity aspect_ratio     angle convexity
#> 10085 10084 238.4863  1.085386    0.7328259     1.469689 139.27116 0.9920712
#> 10176 10175 493.0881  1.067844    0.6938556     1.388669 148.65620 0.9944043
#> 10182 10181 177.5525  1.061164    0.6553132     1.323877 122.36570 0.9914434
#> 10190 10189 323.6471  1.153316    0.8276987     1.781982  45.38620 0.9869685
#> 10247 10246 416.4098  1.072617    0.7292323     1.461420  97.54036 0.9953556
#> 10317 10316 269.0076  1.067543    0.7193451     1.439568  11.75157 0.9958743
#> 10427 10426 375.6754  1.028084    0.1588077     1.012854  16.43842 0.9945204
#> 10433 10432 225.9670  1.038551    0.4909133     1.147831 162.73535 0.9912810
#> 10441 10440 310.9183  1.052642    0.6360212     1.295888  72.76263 0.9936670
#> 10487 10486 329.1423  1.084983    0.7364244     1.478153 159.36537 0.9916191
#>                 geometry sample_id
#> 10085 POINT (5768 14250)     Vis5A
#> 10176 POINT (6909 14910)     Vis5A
#> 10182 POINT (6408 14760)     Vis5A
#> 10190 POINT (6387 13995)     Vis5A
#> 10247 POINT (5830 14742)     Vis5A
#> 10317 POINT (5408 13008)     Vis5A
#> 10427 POINT (6004 14542)     Vis5A
#> 10433 POINT (6299 14453)     Vis5A
#> 10441 POINT (6137 13068)     Vis5A
#> 10487 POINT (6081 13626)     Vis5A
#> 

# Set all annotation geometries, in a named list
toy <- readRDS(system.file("extdata/sfe_toy.rds",
    package = "SpatialFeatureExperiment"
))
ag <- readRDS(system.file("extdata/ag.rds",
    package = "SpatialFeatureExperiment"
))
annotGeometries(toy) <- list(hull = ag)

# Get names of annotation geometries
annotGeometryNames(sfe_small)
#> [1] "tissueBoundary"      "myofiber_full"       "myofiber_simplified"
#> [4] "nuclei"              "nuclei_centroid"    

# Set names of annotation geometries
annotGeometryNames(toy) <- "foo"

# Get a specific annotation geometry by name
# sample_id is optional when there is only one sample present
nuclei <- annotGeometry(sfe_small, type = "nuclei", sample_id = "Vis5A")

# Get a specific annotation geometry by index
tb <- annotGeometry(sfe_small, type = 1L)

# Set a specific annotation geometry
annotGeometry(sfe_small, type = "nuclei2") <- nuclei

# Special convenience function for tissue boundaries
# Getter
tb <- tissueBoundary(sfe_small, sample_id = "Vis5A")
# Setter
tissueBoundary(sfe_small, sample_id = "Vis5A") <- tb
```
