#' SpatialFeatureExperiment coercion methods
#'
#' The \code{SpatialFeatureExperiment} class inherits from
#' \code{SpatialExperiment}, which in turn inherits from
#' \code{SingleCellExperiment}. A \code{SpatialExperiment} object with
#' geometries in \code{colGeometries} in the \code{int_colData},
#' \code{rowGeometries} in the \code{int_elementMetadata}, or
#' \code{annotGeometries} in the \code{int_metadata} can be directly converted
#' to \code{SpatialFeatureExperiment} with \code{as(spe,
#' "SpatialFeatureExperiment")}. A \code{SpatialExperiment} object without the
#' geometries can also be converted; the coordinates in the \code{spatialCoords}
#' field will be used to make POINT geometries named "centroids" to add to
#' \code{colGeometries}. The geometries can also be supplied separately when
#' using \code{toSpatialFeatureExperiment}. Images are converted to \code{SpatRaster}.
#'
#' @inheritParams SpatialFeatureExperiment
#' @inheritParams SpatialExperiment::toSpatialExperiment
#' @inheritParams readVizgen
#' @param x A \code{SpatialExperiment} or \code{Seurat} object to be coerced to a
#'   \code{SpatialFeatureExperiment} object.
#' @param image_scalefactors # A \code{character}, choose between "lowres" or "hires".
#'  Only for 10X Visium, image scaling factors are from `scalefactors_json.json` file.
#' @param unit # Default unit is \code{"micron"}. However for Visium one can choose 
#'  between \code{"micron"} or \code{"full_res_image_pixel"}.
#' @param BPPARAM Deprecated when coercing from \code{SpatialExperiment},
#'  but is used when coercing from \code{Seurat} object.
#' @return A \code{SpatialFeatureExperiment} object
#' @importFrom S4Vectors make_zero_col_DFrame
#' @importFrom SpatialExperiment spatialCoords toSpatialExperiment
#' @name SpatialFeatureExperiment-coercion
#' @aliases toSpatialFeatureExperiment
#' @concept SpatialFeatureExperiment class
#' @examples
#' library(SpatialExperiment)
#' example(read10xVisium)
#' # There can't be duplicate barcodes
#' colnames(spe) <- make.unique(colnames(spe), sep = "-")
#' rownames(spatialCoords(spe)) <- colnames(spe)
#' sfe <- toSpatialFeatureExperiment(spe)
#' # For coercing Seurat to SFE see this -> ./vignettes/seurat_sfe_coerce.Rmd
NULL

.sc2cg <- function(coords_use, spotDiameter = NA) {
    if (is.null(colnames(coords_use)))
        colnames(coords_use) <- paste0("V", seq_len(ncol(coords_use)))
    cg_sfc <- df2sf(coords_use, spatialCoordsNames = colnames(coords_use),
                    spotDiameter = spotDiameter,
                    geometryType = "POINT")
    rownames(cg_sfc) <- rownames(coords_use)
    cg_sfc
}
setAs(
    from = "SpatialExperiment", to = "SpatialFeatureExperiment",
    function(from) {
      cg <- int_colData(from)[["colGeometries"]]
        if (is.null(cg)) {
            coords_use <- spatialCoords(from)
            if (is.null(rownames(coords_use))) {
                rownames(coords_use) <- colnames(from)
            }
            cg <- .sc2cg(coords_use)
            int_colData(from)[["colGeometries"]] <-
                make_zero_col_DFrame(nrow(int_colData(from)))
            int_colData(from)$colGeometries$centroids <- cg
            from
        }
        .spe_to_sfe(from, int_colData(from)[["colGeometries"]],
            int_elementMetadata(from)[["rowGeometries"]],
            int_metadata(from)[["annotGeometries"]],
            spatialCoordsNames(from), "POLYGON",
            int_metadata(from)[["spatialGraphs"]],
            spotDiameter = NA,
            int_metadata(from)[["unit"]]
        )
    }
)

setAs(from = "SingleCellExperiment", to = "SpatialFeatureExperiment",
      function(from) {
          spe <- as(from, "SpatialExperiment")
          as(spe, "SpatialFeatureExperiment")
      })

#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod(
    "toSpatialFeatureExperiment", "SpatialExperiment",
    function(x, colGeometries = NULL, rowGeometries = NULL,
             annotGeometries = NULL, spatialCoordsNames = c("x", "y"),
             annotGeometryType = "POLYGON",
             spatialGraphs = NULL, spotDiameter = NA, unit = NULL,
             BPPARAM = deprecated()) {
        if (is_present(BPPARAM)) {
            deprecate_warn("1.6.0", "toSpatialFeatureExperiment(BPPARAM = )")
        }
        if (is.null(colGeometries)) {
            colGeometries <- int_colData(x)$colGeometries
        }
        if (is.null(rowGeometries)) {
            rowGeometries <- int_elementMetadata(x)$rowGeometries
        }
        if (is.null(annotGeometries)) {
            annotGeometries <- int_metadata(x)$annotGeometries
        }
        if (is.null(spatialGraphs)) {
            spatialGraphs <- int_metadata(x)$spatialGraphs
        }
        .spe_to_sfe(
            x, colGeometries, rowGeometries, annotGeometries,
            spatialCoordsNames, annotGeometryType,
            spatialGraphs, spotDiameter, unit
        )
    }
)

#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod("toSpatialFeatureExperiment", "SingleCellExperiment",
          function(x, sample_id="sample01",
                   spatialCoordsNames = c("x", "y"),
                   spatialCoords=NULL,
                   colGeometries = NULL, rowGeometries = NULL,
                   annotGeometries = NULL,
                   annotGeometryType = "POLYGON",
                   spatialGraphs = NULL, spotDiameter = NA,
                   scaleFactors=1,
                   imageSources=NULL,
                   image_id=NULL,
                   loadImage=TRUE,
                   imgData=NULL,
                   unit = NULL,
                   BPPARAM = deprecated()) {
              if (is_present(BPPARAM)) {
                  deprecate_warn("1.6.0", "toSpatialFeatureExperiment(BPPARAM = )")
              }
              spe <- toSpatialExperiment(x, sample_id=sample_id,
                                         spatialCoordsNames=spatialCoordsNames,
                                         spatialCoords=spatialCoords,
                                         scaleFactors=scaleFactors,
                                         imageSources=imageSources,
                                         image_id=image_id,
                                         loadImage=loadImage,
                                         imgData=imgData)
              toSpatialFeatureExperiment(spe, colGeometries = colGeometries,
                                         rowGeometries = rowGeometries,
                                         annotGeometries = annotGeometries,
                                         spatialCoordsNames = spatialCoordsNames,
                                         annotGeometryType = annotGeometryType,
                                         spatialGraphs = spatialGraphs,
                                         spotDiameter = spotDiameter,
                                         unit = unit)
          })

# Converter function from Seurat (v4 & v5) object to SpatialFeatureExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom methods slot
#' @importFrom SingleCellExperiment SingleCellExperiment mainExpName altExp altExpNames reducedDim
.seu_to_sfe <- 
  function(seu_obj = NULL,
           add_molecules = TRUE,
           flip = c("geometry", "image", "none"),
           image_scalefactors = c("lowres", "hires"),
           unit = NULL,
           BPPARAM = SerialParam()) 
  { 
    # issue message for packages that need to be installed a priori
    check_installed(c("tidyverse", "sf", "BiocParallel", 
                      "Matrix", "Seurat", "SeuratObject", "sfheaders"))
    # checks which Seurat version is present ----
    seu_version <- 
      Biobase::package.version("Seurat") |> 
      grep("5", x =_, value = TRUE)
    
    # support LayerData if Seurat is v5
    .GetCounts <-
      ifelse(length(seu_version) != 0, LayerData, GetAssayData)
    # internal metadata getter 
    .getMeta <- function(object = NULL) {
      if (is(object, "Seurat")) {
        return(slot(object, name = "meta.data"))
      } else {
        return(colData(object) |> 
                 slot(name = "listData") |> 
                 as.data.frame.list())
      }
    }
    # set internal getter for Seurat obj FOV
    .Images <- SeuratObject::Images
    
    # Convert from Seurat to SFE ===== ####
    # TODO (optional) support non-spatial Seurat to SFE as well? ----
    
    if (!is.null(seu_obj)) {
      flip <- match.arg(flip)
      # use existing Assays
      assays_name <- SeuratObject::Assays(seu_obj)
      
      # TODO: add support when Seurat has multiple FOVs/tissue sections ----
      # check FOVs
      if (.Images(seu_obj) |> length() > 1) {
        message(">>> ", .Images(seu_obj) |> length(), 
                " FOVs are found, each will be used as `sample_id`: ", 
                "\n", paste0(.Images(seu_obj), "\n"))
      }
      
      # loop for multiple FOVs/tissue sections ----
      obj_list <-
        lapply(.Images(seu_obj) |> seq(), function(fov_section) {
          message(">>> Seurat Assays found: ", paste0(assays_name, collapse = ", "), "\n",
                  ">>> ", DefaultAssay(seu_obj), " -> will be used as 'Main Experiment'")
          # Make `sf` df geometries
          # get FOVs names
          fovs <- seu_obj[[.Images(seu_obj)[fov_section]]] |> names()
          if(!is.null(fovs)) {
            geoms <- 
              match.arg(fovs, c("centroids", "segmentation", "molecules", "boxes"),
                        several.ok = TRUE)
          } else { geoms <- NULL }
          
          # Support if Seurat object is Visium
          is_Visium <- grep("Visium", 
                            class(seu_obj[[.Images(seu_obj)[fov_section]]]), value = TRUE)
          if (!length(is_Visium) == 0) {
            # get Visium spots coords
            message(">>> Seurat spatial object found: ", is_Visium, "\n",
                    ">>> 'full_res_image_pixel' units will be used ->", "\n",
                    "ie 'imagerow' & 'imagecol' without scaling factors", "\n",
                    ">>> set `unit = 'micron'` to convert spot coordinates to micron space")
            # get image/spots info
            meta_df_image <- 
              slot(seu_obj, name = "images")[[fov_section]] |> 
              slot(name = "coordinates")
            # get xy coords in "pixel" space
            # NOTE: Seurat "spot_diameter_fullres" from `scalefactors_json.json`
            # -> https://github.com/satijalab/seurat/blob/develop/R/preprocessing.R#L1175
            spots <-
              # select coords
              meta_df_image[,grep("image", names(meta_df_image))] |>
              as.matrix()
            colnames(spots)[grep("row", colnames(spots))] <- "pxl_row_in_fullres"
            colnames(spots)[grep("col", colnames(spots))] <- "pxl_col_in_fullres"
            # keep additional metadata
            meta_df_image <- 
              meta_df_image[,grep("tissue|^row|^col", names(meta_df_image))]
            colnames(meta_df_image)[grep("row", colnames(meta_df_image))] <- "array_row"
            colnames(meta_df_image)[grep("col", colnames(meta_df_image))] <- "array_col"
            meta_df_image$in_tissue <- ifelse(meta_df_image$tissue == 1, TRUE, FALSE)
            meta_df_image$tissue <- NULL
            spot_diameter <-
              slot(seu_obj, name = "images")[[fov_section]] |> 
              slot("scale.factors") |> _$spot
          }
          
          # Stop if no FOVs or @images are found
          if (is.null(geoms) && length(is_Visium) == 0) {
            stop(">>> Spatial data is not found, ", "neither FOVs nor `@images` (for Visium)", "\n", 
                 ">>> Check input argument for Seurat object -> `seu_obj` ")
          }
          message(">>> Seurat spatial FOVs found: ", paste0("\n", geoms), "\n",
                  ">>> Generating `sf` geometries")
          
          # Set cell IDs for each FOV/tissue section
          cell_ids_fov <- 
            SeuratObject::Cells(seu_obj[[.Images(seu_obj)[fov_section]]])
          
          # Get centroids
          if (any(geoms == "centroids")) {
            cents <- 
              seu_obj[[.Images(seu_obj)[fov_section]]][["centroids"]] |> 
              as(object = _, "sf")
            rownames(cents) <- cell_ids_fov
          } else { cents <- NULL }
          
          # Get cell polygon geometries
          if (any(geoms == "segmentation")) {
            polys <- 
              seu_obj[[.Images(seu_obj)[fov_section]]][["segmentation"]]
          } else if (!any(geoms == "segmentation") &&
                     any(geoms == "boxes")) {
            # if no segs, use boxes if present
            polys <- 
              seu_obj[[.Images(seu_obj)[fov_section]]][["boxes"]]
          } else { polys <- NULL }
          if (!is.null(polys)) {
            polys <- as(polys, "sf")
            # add cell ids
            polys$ID <- cell_ids_fov
            rownames(polys) <- cell_ids_fov
          }
          
          # Get molecules if provided
          if (add_molecules && any(geoms == "molecules")) {
            mols <- seu_obj[[.Images(seu_obj)[fov_section]]][["molecules"]]
            mols <-
              bplapply(seq(mols),
                       function(i) {
                         data.table::data.table(ID = names(mols)[i], 
                                    slot(mols[[i]], "coords"))
                         }, BPPARAM = BPPARAM) |> rbindlist()
            mols <- df2sf(mols, geometryType = "MULTIPOINT", group_col = "ID")
            } else { mols <- NULL }
          
          # Get sparse count matrices
          # prepare assays
          assays_sfe <- 
            list(counts = .GetCounts(seu_obj, "counts",
                                     assay = DefaultAssay(seu_obj)),
                 logcounts = .GetCounts(seu_obj, "data",
                                        assay = DefaultAssay(seu_obj)),
                 scaledata = .GetCounts(seu_obj, "scale.data",
                                        assay = DefaultAssay(seu_obj)))
          # remove empty elements, if eg "scale.data" is not present
          assays_sfe <- assays_sfe[lapply(assays_sfe, length) > 0] 
          # # using cell IDs from FOV/tissue section cell_ids_fov
          assays_sfe <- 
            lapply(seq(assays_sfe), function(i) assays_sfe[[i]][,cell_ids_fov]) |>
            # keep names
            setNames(object = _, names(assays_sfe))
          
          assays_sfe_out <<- assays_sfe
          
          # Make SFE from DefaultAssay of Seurat ----
          message("\n", ">>> Creating `SFE` object -> ", .Images(seu_obj)[fov_section])
          sfe <- SpatialFeatureExperiment(assays = assays_sfe,
                                          colData = 
                                            if (!length(is_Visium) == 0) {
                                              cbind(.getMeta(seu_obj)[cell_ids_fov,],
                                                    meta_df_image)
                                            } else { .getMeta(seu_obj)[cell_ids_fov,] },
                                          sample_id = gsub("_|-|[.]", "_", 
                                                           x = .Images(seu_obj)[fov_section]),
                                          colGeometries = 
                                            if (!is.null(cents))
                                              # use centroids POINT here
                                              list(centroids = cents) 
                                            else if (!length(is_Visium) == 0)
                                              # spot centroids for Visium
                                              list("spotPoly" = 
                                                     .sc2cg(spots, 
                                                            spotDiameter = spot_diameter)),
                                          unit = 
                                            if (!is.null(unit) &&
                                                is.character(unit)) 
                                              unit
                                            #..for Visium only
                                            else if (is.null(unit) && 
                                                     !length(is_Visium) == 0)
                                              "full_res_image_pixel" else "micron",
                                          mainExpName = DefaultAssay(seu_obj))
          # add sample_id to centroids
          colGeometry(sfe, 1)$sample_id <- sampleIDs(sfe)
          
          sfe_out <<- sfe
          
          # TODO: add reducedDim ----
          dimRed <- slot(seu_obj, "reductions") |> names()
          message(">>> Adding Dimensionality Reduction: ", paste0(dimRed, collapse = ", "))
          for (i in seq(dimRed))
            reducedDim(x = sfe, 
                       type = toupper(dimRed[i])) <- SeuratObject::Embeddings(sfe, reduction = dimRed[i])
          # add variance and other parts?
          

          # flip geometry both for col & row geoms ----
          if (flip == "geometry") {
            # flip the coordinates
            mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
            if (!is.null(polys))
                st_geometry(polys) <- st_geometry(polys) * mat_flip
            if (!is.null(mols))
              st_geometry(mols) <- st_geometry(mols) * mat_flip
            if (!is.null(cents))
              st_geometry(colGeometry(sfe , "centroids")) <-
                st_geometry(colGeometry(sfe , "centroids")) * mat_flip
            }
          
          # add polygon geometries ----
          if (!is.null(polys)) {
            # sanity on geometries
            polys <- .check_st_valid(polys)
            cellSeg(sfe) <- polys
            cellSeg(sfe)$sample_id <- sampleIDs(sfe)
          }
          
          # convert pixels to microns (Visium only) ----
          if (unit(sfe) == "micron" && !length(is_Visium) == 0) {
            message(">>> Converting pixels to microns")
            scale_fct <- .pixel2micron(sfe)
            cg <- spotPoly(sfe)
            cg$geometry <- cg$geometry * scale_fct
            spotPoly(sfe) <- cg
            spotPoly(sfe)$sample_id <- sampleIDs(sfe)
            
            # TODO (under construction): add Images to SFE from Seurat or from custom path ----
            # optionally make arg to provide image path to add image, ie for Xenium, Vizgen &| Visium
            if (FALSE) {
              # if Seurat is Visium
              # Scale factors for images
              # note, scale_imgs is a value in, eg @scale.factors$lowres
              scale_imgs <- 
                slot(seu_obj, name = "images")[[fov_section]] |> 
                slot(name = "scale.factors") |> _[image_scalefactors]
              scale_imgs <- scale_imgs / scale_fct
              
              # eg image part
              slot(seu_obj, name = "images")[[fov_section]] |> 
                slot(name = "image") |> str()
              # convert to image raster terra
              
              # add to sfe
              
              # eg
              # Set up ImgData
              img_fns2 <- file.path(dirs[i], "spatial", img_fns)
              img_dfs <- lapply(seq_along(img_fns), function(j) {
                .get_imgData(img_fns2[j], sample_id = sample_id[i],
                             image_id = names(img_fns)[j],
                             extent = NULL, scale_fct = scale_imgs[[j]],
                             flip = TRUE)
              })
              img_df <- do.call(rbind, img_dfs)
              imgData(o) <- img_df
            }
            
          }
          
          # add transcript coordinates ----
          if (add_molecules && !is.null(mols)) {
            # store molecules in rowGeometries
            message(">>> Storing molecule coordinates in `rowGeometries`")
            # sanity, if number of genes doesn’t overlap, subset SFE to genes present in mols
            genes_mols <- unique(mols$ID)
            if (!all(genes_mols %in% rownames(sfe))) {
              gene_indx <-
                which(rownames(sfe) %in% genes_mols)
              # genes/features that are removed
              genes_rm <- rownames(sfe)[-gene_indx]
              warning(">>> Total number of genes in SFE object does NOT overlap with `transcript` coordinates genes", "\n",
                      ">>> Subsetting object and `transcript` coordinates to keep only matching genes", "Total of ", length(genes_rm),
                      " genes are removed")
              # subset sfe
              sfe <- sfe[gene_indx,]
              # subset transcripts keep on
              mols <- mols[mols$ID %in% rownames(sfe),]
            }
            if (is(mols, "sf")) {
              rownames(mols) <- unique(mols$ID)
              txSpots(sfe, withDimnames = TRUE) <- mols
              # add sample id
              txSpots(sfe)$sample_id <- sampleIDs(sfe)
            }
          }
          
          # TODO (alternatively): use `MultiAssayExperiment` instead? ----
          # Currently using `altExp` if > 1 Seurat Assays are present ----
          if (length(assays_name) > 1) {
            # keep assays minus default one
            assays_name <- setdiff(assays_name, DefaultAssay(seu_obj))
            message(">>> Adding Seurat Assay(s) as Alternative Experiment(s): ", 
                    paste0("\n", assays_name))
            # loop and add multiple assays
            for (a in seq(assays_name)) {
              # prepare assays
              assays_add <- 
                list(counts = .GetCounts(seu_obj, "counts",
                                         assay = assays_name[a]),
                     logcounts = .GetCounts(seu_obj, "data",
                                            assay = assays_name[a]),
                     scaledata = .GetCounts(seu_obj, "scale.data",
                                            assay = assays_name[a]))
              # remove empty elements, if eg "scale.data" 
              assays_add <- assays_add[lapply(assays_add, length) > 0]
              
              # subset to keep cells only for selected FOV/tissue section
              assays_add <- 
                lapply(seq(assays_add), function(i) assays_add[[i]][,cell_ids_fov]) |>
                # keep names
                setNames(object = _, names(assays_add))
              altExp(sfe, assays_name[a]) <- 
                SingleCellExperiment(assays = assays_add) |> 
                toSpatialExperiment(sample_id = gsub("_|-|[.]", "", x = .Images(seu_obj)[fov_section])) |>
                toSpatialFeatureExperiment(x = _,
                                           spatialCoordsNames = c("x", "y"), 
                                           colGeometries = list(centroids = colGeometry(sfe)), 
                                           unit = unit(sfe))
            }
          }
          return(sfe)
        })
      
      obj_list_test <<- obj_list
      
      if (length(obj_list) > 1) {
        message(">>> Combining ", length(obj_list), 
                " SFE object(s) with unique `sample_id`")
        
        # TODO: cbind SFEs and altExps ----
        # extract altExp list -> lapply, add sample_id to each altExp
        # cbind SFEs without altExps
        # cbind altExps
        # add combined altExps to combined SFEs
        
        
        
        # TODO: return single combined object (when multiple samples/FOVs) ----
        #return(obj_list)
      } else if (length(obj_list) == 1){
        message("\n", ">>> `SFE` object is ready!")
        return(obj_list[[1]])
      }
    } else { stop("No Seurat object in `x` was found") }
  }

# Set formal method for Seurat to SFE
#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod("toSpatialFeatureExperiment", "Seurat",
          definition = 
            function(x,
                     add_molecules = TRUE,
                     flip = c("geometry", "image", "none"),
                     image_scalefactors = c("lowres", "hires"),
                     unit = NULL,
                     BPPARAM = SerialParam()) {
              .seu_to_sfe(seu_obj = x,
                          add_molecules = add_molecules,
                          flip = flip,
                          image_scalefactors = image_scalefactors,
                          unit = unit,
                          BPPARAM = BPPARAM)
            })
