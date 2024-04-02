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
#' @param x A \code{SpatialExperiment} object to be coerced to a
#'   \code{SpatialFeatureExperiment} object.
#' @param BPPARAM Deprecated.
#' @return An SFE object
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
             spatialGraphs = NULL, spotDiameter = NA, unit = NULL) {
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
                   unit = NULL) {
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


#' Coerce from Seurat to SpatialFeatureExperiment
#' This function converts \code{Seurat} (v4 & v5) object to \code{SpatialFeatureExperiment} object
#' @inheritParams readVizgen
#' @param seurat_obj A \code{Seurat} object.
#' @param image_scalefactors # A \code{character}, choose between "lowres" or "hires".
#'  Only for 10X Visium, image scaling factors are from `scalefactors_json.json`.
#' @param unit_use # Default unit is \code{"micron"}. However for Visium one can choose 
#'  between \code{"micron"} or \code{"full_res_image_pixel"}.
#' @return A \code{SpatialFeatureExperiment} object
#' @export
#' @importFrom SummarizedExperiment assays
#' @importFrom methods slot
#' @importFrom SingleCellExperiment SingleCellExperiment mainExpName altExp altExpNames
#' @importFrom BiocParallel bplapply
#' 
#' @examples 
#' # TODO or see ./vignettes/seurat_sfe_coerce.Rmd
#' 

# TODO, change dplyr functions, eg mutate, to base ----

.seurat_to_sfe <- 
  function(seurat_obj = NULL,
           add_molecules = TRUE,
           flip = c("geometry", "image", "none"),
           image_scalefactors = "lowres",
           unit_use = NULL,
           BPPARAM = SerialParam()) 
  { # issue message for packages that need to be installed a priori
    pkgs <- c("tidyverse", "sf", "sp", "BiocParallel", 
              "Matrix", "Seurat", "SeuratObject")
    pkgs_check <-
      lapply(pkgs |> length() |> seq(), function(i) 
      { !requireNamespace(pkgs[i], quietly = TRUE) }) |> unlist()
    { if (any(which(pkgs_check) > 0)) 
    { message("Please install ->", "\n",
              paste0("'", pkgs[which(pkgs_check)], "'", collapse = ", "), " for this function")}
      }
    
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
    
    #### ===== Convert from Seurat to SFE ===== ####
    # TODO (optional) support non-spatial Seurat to SFE as well? ----
    
    if (!is.null(seurat_obj)) {
      # use existing Assays
      assays_name <- SeuratObject::Assays(seurat_obj)
      
      # TODO: add support when Seurat has multiple FOVs/tissue sections ----
      # loop/lapply
      # seurat_obj[[Images(seurat_obj)[1]]] # coords or FOVs of 1st section
      # implement:
      # extract FOVs
      # subset single cell object (no FOVs) using cell id from "centoirds" fov
      # similar for Visium -> slot(seurat_obj, name = "images")[[1]] or seurat_obj[[Images(seurat_obj)[1]]]
      # add sample_id <- "" another value to cbind the objs, same for for altExp
      # or, return the list() of objects
      
      # check FOVs
      if (Images(seurat_obj) |> length() > 1) {
        message(">>> ", Images(seurat_obj) |> length(), 
                " FOVs are found, each will be used as `sample_id`: ", 
                "\n", paste0(Images(seurat_obj), "\n"))
      }
      
      # loop for multiple FOVs/tissue sections ----
      obj_list <-
        lapply(Images(seurat_obj) |> seq(), function(wooh) {
          message(">>> Seurat Assays found: ", paste0(assays_name, ", "), "\n",
                  ">>> ", DefaultAssay(seurat_obj), " -> will be used as 'Main Experiment'")
          # Make `sf` df geometries
          # get FOVs names
          fovs <- seurat_obj[[Images(seurat_obj)[wooh]]] |> names()
          if(!is.null(fovs)) {
            geoms <- 
              match.arg(fovs, c("centroids", "segmentation", "molecules", "boxes"),
                        several.ok = TRUE)
          } else { geoms <- NULL }
          
          # Support if Seurat object is Visium
          is_Visium <- grep("Visium", 
                            class(seurat_obj[[Images(seurat_obj)[wooh]]]), value = TRUE)
          if (!length(is_Visium) == 0) {
            # get Visium spots coords
            message(">>> Seurat spatial object found: ", is_Visium, "\n",
                    ">>> 'full_res_image_pixel' units will be used ->", "\n",
                    "ie 'imagerow' & 'imagecol' without scaling factors", "\n",
                    ">>> set `unit_use = 'micron'` to convert spot coordinates to micron space")
            # get image/spots info
            meta_df_image <- 
              slot(seurat_obj, name = "images")[[wooh]] |> 
              slot(name = "coordinates")
            # get xy coords in "pixel" space
            # NOTE: Seurat "spot_diameter_fullres" from `scalefactors_json.json`
            # -> https://github.com/satijalab/seurat/blob/develop/R/preprocessing.R#L1175
            spots <-
              meta_df_image |>
              # select coords
              dplyr::transmute(pxl_col_in_fullres = imagecol,
                               pxl_row_in_fullres = imagerow) |>
              as.matrix()
            # keep additional metadata
            meta_df_image <- 
              meta_df_image |>
              dplyr::transmute(in_tissue = case_when(tissue == 1 ~ TRUE, 
                                                     tissue != 1 ~ FALSE),
                               array_row = row, 
                               array_col = col)
            spot_diameter <-
              slot(seurat_obj, name = "images")[[wooh]] |> 
              slot("scale.factors") |> _$spot
          }
          
          # Stop if no FOVs or @images are found
          if (is.null(geoms) && length(is_Visium) == 0) {
            stop(">>> Spatial data is not found, ", "neither FOVs nor `@images` (for Visium)", "\n", 
                 ">>> Check input argument for Seurat object -> `seurat_obj` ")
          }
          message(">>> Seurat spatial FOVs found: ", paste0("\n", geoms), "\n",
                  ">>> Generating `sf` geometries")
          
          # Set cell IDs for each FOV/tissue section
          cell_ids_fov <- seurat_obj[[Images(seurat_obj)[wooh]]] |> Cells()
          
          # Get centroids
          if (any(geoms == "centroids")) {
            cents <- 
              seurat_obj[[Images(seurat_obj)[wooh]]][[grep("centroids", geoms, value = TRUE)]] |> 
              as(object = _, "sf")
            rownames(cents) <- cell_ids_fov
          } else { cents <- NULL }
          
          # Get cell polygon geometries
          if (any(geoms == "segmentation")) {
            polys <- 
              seurat_obj[[Images(seurat_obj)[wooh]]][[grep("segmentation", geoms, value = TRUE)]] |>
              slot(name = "polygons")
          } else if (!any(geoms == "segmentation") &&
                     any(geoms == "boxes")) {
            # if no segs, use boxes if present
            polys <- 
              seurat_obj[[Images(seurat_obj)[wooh]]][[grep("boxes", geoms, value = TRUE)]] |>
              slot(name = "polygons")
          } else { polys <- NULL }
          if (!is.null(polys)) {
            polys <- 
              bplapply(polys |> seq(), function(i) {
                st_sf(ID = names(polys)[i],
                      geometry = polys[[i]] |> 
                        list() |>
                        sp::SpatialPolygons(Srl = _) |> 
                        as(object = _, "sf") |> st_geometry())
              }, BPPARAM = BPPARAM) |>
              do.call(bind_rows, args = _)
            rownames(polys) <- polys$ID
          }
          
          # Get molecules if provided
          if (add_molecules && any(geoms == "molecules")) {
            mols <- seurat_obj[[Images(seurat_obj)[wooh]]][[grep("molecules", geoms, value = TRUE)]]
            mols <- 
              lapply(mols |> seq(), function(i) {
                st_sf(geometry = mols[[i]] |>
                        list() |>
                        sp::SpatialMultiPoints(coords = _) |>
                        as(object = _, "sf")) |>
                  mutate(ID = names(mols)[i])}) |> 
              do.call(bind_rows, args = _)
          } else { mols <- NULL }
          
          # Get sparse count matrices
          # prepare assays
          assays_sfe <- 
            list(counts = .GetCounts(seurat_obj, "counts",
                                     assay = DefaultAssay(seurat_obj)),
                 logcounts = .GetCounts(seurat_obj, "data",
                                        assay = DefaultAssay(seurat_obj)),
                 scaledata = .GetCounts(seurat_obj, "scale.data",
                                        assay = DefaultAssay(seurat_obj)))
          # remove empty elements, if eg "scale.data" is not present
          # alternative pass it |> purrr::compact(.)
          assays_sfe <- assays_sfe[lapply(assays_sfe, length) > 0] 
          # # using cell IDs from FOV/tissue section cell_ids_fov
          assays_sfe <- 
            lapply(seq(assays_sfe), function(i) assays_sfe[[i]][,cell_ids_fov]) |>
            # keep names
            setNames(object = _, names(assays_sfe))
          
          assays_sfe_out <<- assays_sfe
          
          # Make SFE from DefaultAssay of Seurat ----
          message(">>> Creating `SFE` object -> ", Images(seurat_obj)[wooh])
          sfe <- SpatialFeatureExperiment(assays = assays_sfe,
                                          colData = 
                                            if (!length(is_Visium) == 0) {
                                              cbind(.getMeta(seurat_obj)[cell_ids_fov,],
                                                    meta_df_image)
                                            } else { .getMeta(seurat_obj)[cell_ids_fov,] },
                                          sample_id = gsub("_|-|[.]", "_", 
                                                           x = Images(seurat_obj)[wooh]),
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
                                            if (!is.null(unit_use) &&
                                                is.character(unit_use)) 
                                              unit_use
                                          #..for Visium only
                                          else if (is.null(unit_use) &&
                                                   !length(is_Visium) == 0)
                                            "full_res_image_pixel" else "micron",
                                          mainExpName = DefaultAssay(seurat_obj))
          # add sample_id to centroids
          colGeometry(sfe, 1)$sample_id <- sampleIDs(sfe)
          
          sfe_out <<- sfe
          polys_out <<- polys
          mols_out <<- mols
          if (!length(is_Visium) == 0)
            spots_out <<- spots
          
          # add polygon geometries ----
          if (!is.null(polys)) {
            # sanity on geometries
            polys <- .check_st_valid(polys)
            # flip geometry
            if (flip == "geometry" && !is.null(polys)) {
              # Flip the coordinates
              mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
              st_geometry(polys) <- st_geometry(polys) * mat_flip
            }
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
            
            # TODO: add Images to SFE from Seurat or custom path ----
            # for now if Seurat is Visium
            # optionally image path could be provided and image added, eg for Xenium, Vizgen
            #..under construction
            if (FALSE) {
              
              # Scale factors for images
              # TODO: scale_imgs is a value in, eg @scale.factors$lowres ----
              scale_imgs <- 
                slot(seurat_obj, name = "images")[[wooh]] |> 
                slot(name = "scale.factors") |> _[image_scalefactors]
              scale_imgs <- scale_imgs / scale_fct
              
              # eg image part
              slot(seurat_obj, name = "images")[[wooh]] |> 
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
              mols <- mols_out[mols$ID %in% rownames(sfe),]
            }
            if (is(mols, "sf")) {
              rownames(mols) <- unique(mols$ID)
              mols_out <<- mols
              txSpots(sfe, withDimnames = TRUE) <- mols
              # add sample id
              txSpots(sfe)$sample_id <- sampleIDs(sfe)
            }
          }
          
          # TODO: if multiple Assays are present ----
          # Add as Alternative Experiment
          # unless it's Xenium data with these assays:
          # 'Gene Expression''Negative Control Probe''Negative Control Codeword''Unassigned Codeword'
          # then all Assays gene will go into 1 assay with all the genes
          
          # eg --->>>  
          # rowData(sfe, 2) |> str
          # Formal class 'DFrame' [package "S4Vectors"] with 6 slots
          #..@ rownames       : chr [1:541] "ABCC9" "ADAM17" "ADAMTS12" "ADAMTS16" ...
          #..@ nrows          : int 541
          #..@ elementType    : chr "ANY"
          #..@ elementMetadata: NULL
          #..@ metadata       : list()
          #..@ listData       :List of 3
          #.. ..$ ID    : chr [1:541] "ENSG00000069431" "ENSG00000151694" "ENSG00000151388" "ENSG00000145536" ...
          #.. ..$ Symbol: chr [1:541] "ABCC9" "ADAM17" "ADAMTS12" "ADAMTS16" ...
          #.. ..$ Type  : chr [1:541] "Gene Expression" "Gene Expression" "Gene Expression" "Gene Expression" ...
          # eg <<<---
          
          # Currently adding as Alternative Experiment
          if (length(assays_name) > 1) {
            # keep assays minus default one
            assays_name <- setdiff(assays_name, DefaultAssay(seurat_obj))
            message(">>> Adding Seurat Assay(s) as Alternative Experiment(s): ", 
                    paste0("\n", assays_name))
            # loop and add multiple assays
            for (a in seq(assays_name)) {
              # prepare assays
              assays_add <- 
                list(counts = .GetCounts(seurat_obj, "counts",
                                         assay = assays_name[a]),
                     logcounts = .GetCounts(seurat_obj, "data",
                                            assay = assays_name[a]),
                     scaledata = .GetCounts(seurat_obj, "scale.data",
                                            assay = assays_name[a]))
              # remove empty elements, if eg "scale.data" 
              assays_add <- assays_add[lapply(assays_add, length) > 0]
              
              # subset to keep cells only for selected FOV/tissue section
              assays_add <- 
                lapply(seq(assays_add), function(i) assays_add[[i]][,cell_ids_fov]) |>
                # keep names
                setNames(object = _, names(assays_add))
              altExp(sfe, assays_name[a]) <- 
                SingleCellExperiment::SingleCellExperiment(assays = assays_add) |> 
                SpatialExperiment::toSpatialExperiment(sample_id = gsub("_|-|[.]", "", 
                                                                        x = Images(seurat_obj)[wooh])) |>
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
        
        # TODO: cbind SFEs and altExp ----
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
    }
  }

#' Set formal method for Seurat to SFE
#' @rdname SpatialFeatureExperiment-coercion
#' @export
setMethod("toSpatialFeatureExperiment", "Seurat",
          function(x,
                   add_molecules = TRUE,
                   flip = c("geometry", "image", "none"),
                   image_scalefactors = "lowres",
                   unit_use = NULL,
                   BPPARAM = SerialParam()) {
            .seurat_to_sfe(seurat_obj = x,
                           add_molecules = add_molecules,
                           flip = flip,
                           image_scalefactors = image_scalefactors,
                           unit_use = unit_use,
                           BPPARAM = BPPARAM)
            })

