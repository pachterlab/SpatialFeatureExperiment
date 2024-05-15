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
#' @importFrom methods slot<-
#' @importFrom SingleCellExperiment altExp<- reducedDim<-
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

.sc2cg <- function(coords_use, spotDiameter = NA, ...) {
    if (is.null(colnames(coords_use)))
        colnames(coords_use) <- paste0("V", seq_len(ncol(coords_use)))
    cg_sfc <- df2sf(coords_use, spatialCoordsNames = colnames(coords_use),
                    spotDiameter = spotDiameter,
                    geometryType = "POINT",
                    ... # this arg goes to `sf::st_buffer`, see `df2sf`
                    )
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

#..to get slot or layer names of a selected Assay
.GetSlotNames <- function(object_seu, assay_seu, fov_number) {
    slot_n <-
        Seurat::GetAssay(object_seu, assay_seu)
    if (is(slot_n, "Assay")) {
        # Seurat v4 based object
        slot_n <- slotNames(x = slot_n)[1:2]
    } else if (is(slot_n, "Assay5")) {
        # Seurat v5 based object
        slot_n <- slot(slot_n, name = slotNames(x = slot_n)[1]) |> names()
    }
    # sanity, check slot names, for cases like:
    #c("counts.1", "counts.2", "data.1", "data.2")
    ok_slots <- any(slot_n %in% c("counts", "data"))
    if (!ok_slots) {
        if (length(slot_n) > 2) {
            slot_n <-
                grep(as.character(fov_number), slot_n, value = TRUE)
        } else if (length(slot_n) == 2) {
            slot_n <- slot_n[fov_number]
        }
    }
    return(list(slot_n = slot_n,
                slots_ok = ok_slots))
}

# internal metadata getter for Seurat and SFE objects
.getMeta <- function(object = NULL) {
    if (is(object, "Seurat")) {
        return(slot(object, "meta.data"))
    } else {
        return(colData(object) |>
                   as.data.frame())
    }
}

# Converter function from Seurat (v4 & v5) object to SpatialFeatureExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom methods slot slotNames
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
    check_installed(c("dplyr", "tidyr", "Seurat", "SeuratObject"))
    # checks which Seurat version is present
    seu_version <-
      Biobase::package.version("Seurat") |>
      grep("5", x =_, value = TRUE)

    # Internal getters for Seurat object:
    #..to support LayerData if Seurat is v5
    .GetCounts <-
      ifelse(length(seu_version) != 0,
             SeuratObject::LayerData,
             SeuratObject::GetAssayData)

    #..to get default or main Assay
    .DefaultAssay <- SeuratObject::DefaultAssay

    # Convert from Seurat to SFE ===== ####
    # TODO (optional) support non-spatial Seurat to SFE as well? ----

    if (!is.null(seu_obj)) {
      flip <- match.arg(flip)
      # use existing Assays
      assays_name <- SeuratObject::Assays(seu_obj)

      # Supports when Seurat has multiple FOVs/tissue sections ----
      fov_names <- SeuratObject::Images(seu_obj)
      if (length(fov_names) > 1) {
        message(">>> ", length(fov_names),
                " spatial FOVs are found, each will be used as `sample_id`: ",
                "\n", paste0(fov_names, "\n"))
      }

      # loop for multiple FOVs/tissue sections ----
      # TODO: (enhancement) consider bplapply ----
      #..ie, when looping over FOVs with large number of cells and molecules
      obj_list <-
        lapply(seq(fov_names), function(fov_section) {
          message(">>> Seurat Assays found: ", paste0(assays_name, collapse = ", "), "\n",
                  ">>> ", .DefaultAssay(seu_obj), " -> will be used as 'Main Experiment'")
          # Make `sf` df geometries
          # get spatial S4 FOV object
          fov_seu <- seu_obj[[fov_names[fov_section]]]
          # get FOV main names
          fovs <- names(fov_seu)
          if(!is.null(fovs)) {
            geoms <-
              match.arg(fovs, c("centroids", "segmentation", "molecules", "boxes"),
                        several.ok = TRUE)
          } else { geoms <- NULL }

          # Support when Seurat object is Visium
          is_Visium <- grep("Visium",
                            class(fov_seu), value = TRUE)
          # Stop if FOVs or slot "images" are NOT found
          if (is.null(geoms) && length(is_Visium) == 0) {
            stop(">>> Spatial data is not found, ", "neither FOVs nor slot 'images' (for Visium)", "\n",
                 ">>> Check input argument for Seurat object -> `seu_obj` ")
          }
          if (!length(is_Visium) == 0 && is_Visium == "VisiumV1") {
            # get Visium spots coords
            message(">>> Seurat spatial object found: ", is_Visium, "\n",
                    ">>> 'full_res_image_pixel' units will be used ->", "\n",
                    "ie 'imagerow' & 'imagecol' without scaling factors", "\n",
                    ">>> set `unit = 'micron'` to convert spot coordinates to micron space")
            # get image/spots info
            relevant_slots <- grep("coordinates|scale.factors",
                                   slotNames(fov_seu), value = TRUE)
            # sanity, make sure the relevant slots are present
            if (length(relevant_slots) == 2) {
              meta_df_image <- slot(fov_seu, "coordinates")
              # get xy coords in "pixel" space
              # NOTE: Seurat "spot_diameter_fullres" from `scalefactors_json.json`
              # -> https://github.com/satijalab/seurat/blob/develop/R/preprocessing.R#L1175
              spots <-
                # select coords
                meta_df_image[,grep("image", names(meta_df_image))] |>
                as.matrix()
              colnames(spots)[grep("row", colnames(spots))] <- "pxl_row_in_fullres"
              colnames(spots)[grep("col", colnames(spots))] <- "pxl_col_in_fullres"
              spot_diameter <- slot(fov_seu, "scale.factors")[["spot"]]
              spots <- .sc2cg(spots, spotDiameter = spot_diameter, endCapStyle = "ROUND")
              # keep additional metadata
              meta_df_image <-
                meta_df_image[,grep("tissue|^row|^col", names(meta_df_image))]
              colnames(meta_df_image)[grep("row", colnames(meta_df_image))] <- "array_row"
              colnames(meta_df_image)[grep("col", colnames(meta_df_image))] <- "array_col"
              meta_df_image$in_tissue <- ifelse(meta_df_image$tissue == 1, TRUE, FALSE)
              meta_df_image$tissue <- NULL
              } else { stop(">>> Only ",
                            paste0("'", relevant_slots,"'"),
                            " slot is present in VisiumV1, 'coordinates' and 'scale.factors' are both needed")
              }
            }

          if (!is.null(geoms) && length(is_Visium) == 0)
            message(">>> Seurat spatial FOVs found: ", paste0("\n", geoms))

          message(">>> Generating `sf` geometries")
          # Set cell IDs for each FOV/tissue section
          cell_ids_fov <- SeuratObject::Cells(fov_seu)

          # Get centroids, including VisiumHD
          if (any(geoms == "centroids")) {
            if (!length(is_Visium) == 0 && is_Visium == "VisiumV2") {
              message(">>> Seurat spatial object found: ", is_Visium, " -> VisiumHD", "\n")
              # get image/spots info
              relevant_slots <- grep("boundaries|scale.factors",
                                     slotNames(fov_seu), value = TRUE)
              # sanity, make sure the relevant slots are present
              if (length(relevant_slots) == 2) {
                #meta_df_image <- slot(fov_seu, "coordinates")
                # get xy coords in "pixel" space
                spotsHD <- slot(fov_seu[["centroids"]], "coords")
                rownames(spotsHD) <- cell_ids_fov
                # get spot diameter in pixels
                spot_diameter <- slot(fov_seu, "scale.factors")[["spot"]]
                # get bin size, ie 2, 8 or 16 microns
                bin_um <-
                  sapply(c("2","8","16"), function(i)
                    grep(i, slot(fov_seu, "assay")), USE.NAMES = TRUE) |>
                  unlist() |> names() |> as.integer()
                # converting bins to squre polygons
                message(">>> Making POLYGON geometry for bin size: ", bin_um, " \u03BCm with total of ",
                        length(cell_ids_fov), " cells, this can take some time!", "\n")
                spotsHD <- .sc2cg(spotsHD, spotDiameter = spot_diameter, endCapStyle = "SQUARE")
                cents <- NULL
                } else { stop(">>> Only ",
                            paste0("'", relevant_slots,"'"),
                            " slot is present in VisiumV2, 'coordinates' and 'scale.factors' are both needed")
                }
              } else {
              # centroids will be POINT geometry
              cents <- as(object = fov_seu[["centroids"]], "sf")
              rownames(cents) <- cell_ids_fov
              }
            } else { cents <- NULL }

          # Get cell polygon geometries
          if (any(geoms == "segmentation")) {
            polys <-
              fov_seu[["segmentation"]]
          } else if (!any(geoms == "segmentation") &&
                     any(geoms == "boxes")) {
            # if no segs, use boxes if present
            polys <-
              fov_seu[["boxes"]]
          } else { polys <- NULL }
          if (!is.null(polys)) {
            polys <- as(polys, "sf")
            # add cell ids
            polys$ID <- cell_ids_fov
            rownames(polys) <- cell_ids_fov
          }

          # Get molecules if provided
          if (add_molecules && any(geoms == "molecules")) {
            mols <- fov_seu[["molecules"]]
            mols <-
              bplapply(seq(mols),
                       function(i) {
                         data.table::data.table(ID = names(mols)[i],
                                                slot(mols[[i]], "coords"))
                         }, BPPARAM = BPPARAM) |> rbindlist()
            mols <- df2sf(mols, geometryType = "MULTIPOINT", group_col = "ID")
            } else { mols <- NULL }

          # Get sparse count matrices
          assay_master <-
          if (length(assays_name) == length(fov_names))
            assays_name[fov_section]
          else .DefaultAssay(seu_obj)
          # get slot or layer names
          slot_names <- .GetSlotNames(object_seu = seu_obj,
                                      assay_seu = assay_master,
                                      fov_number = fov_section)
          # prepare assays
          # TODO: this code is a bit clunky, make pretty later
          assays_sfe <-
            list(counts = .GetCounts(seu_obj,
                                     cells = cell_ids_fov,
                                     if (!slot_names[["slots_ok"]] &&
                                         length(slot_names[["slot_n"]]) == 1)
                                       slot_names[["slot_n"]]
                                     else if (!slot_names[["slots_ok"]] &&
                                              length(slot_names[["slot_n"]]) == 2)
                                       slot_names[["slot_n"]][1]
                                     else "counts",
                                     assay = assay_master),
                 logcounts = .GetCounts(seu_obj,
                                        cells = cell_ids_fov,
                                        if (!slot_names[["slots_ok"]] && length(slot_names[["slot_n"]]) == 2)
                                          slot_names[["slot_n"]][2]
                                        else "data",
                                        assay = assay_master),
                 scaledata = .GetCounts(seu_obj,
                                        cells = cell_ids_fov,
                                        "scale.data",
                                        assay = assay_master))
          # remove empty elements, if eg "scale.data" is not present
          assays_sfe <- assays_sfe[lapply(assays_sfe, length) > 0]

          # remove this afterwards
          #assays_sfe_out <<- assays_sfe

          # Make SFE from Seurat ----
          message("\n", ">>> Creating `SFE` object -> ", fov_names[fov_section])
          sfe <- SpatialFeatureExperiment(assays = assays_sfe,
                                          colData =
                                            if (!length(is_Visium) == 0 && is_Visium == "VisiumV1")
                                              cbind(.getMeta(seu_obj)[cell_ids_fov,],
                                                    meta_df_image)
                                            else .getMeta(seu_obj)[cell_ids_fov,],
                                          sample_id = gsub("_|-|[.]", "_",
                                                           x = fov_names[fov_section]),
                                          colGeometries =
                                            if (!length(is_Visium) == 0) {
                                              if (is_Visium == "VisiumV2")
                                                # square bin POLYGON for VisiumHD
                                                list(spotPoly = spotsHD)
                                              else if (is_Visium == "VisiumV1")
                                                # circle spot POLYGON for Visium
                                                list(spotPoly = spots)
                                              } else if (!is.null(cents)) {
                                                # centroids POINT
                                                list(centroids = cents)},
                                          unit =
                                            if (!is.null(unit) && is.character(unit))
                                              unit
                                            # for Visium
                                            else if (is.null(unit) && !length(is_Visium) == 0)
                                              "full_res_image_pixel" else "micron",
                                          mainExpName = assay_master
                                          )
          # add sample_id to centroids
          colGeometry(sfe, 1)$sample_id <- sampleIDs(sfe)

          # remove this afterwards
          #sfe_out <<- sfe

          # flip geometry both for col & row geoms ----
          if (flip == "geometry") {
            # flip the coordinates
            mat_flip <- matrix(c(1,0,0,-1), ncol = 2)
            if (!is.null(polys))
                st_geometry(polys) <- st_geometry(polys) * mat_flip
            if (!is.null(mols))
              st_geometry(mols) <- st_geometry(mols) * mat_flip
            if (!is.null(cents) && length(is_Visium) == 0)
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

          # convert pixels to microns (Visium and VisiumHD only) ----
          if (unit(sfe) == "micron" && !length(is_Visium) == 0) {
            message(">>> Converting pixels to microns")
            if (is_Visium == "VisiumV2")
              # set scaling factor -> microns per pixel
              scale_fct <- bin_um / spot_diameter
            else
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
                slot(seu_obj, "images")[[fov_section]] |>
                slot("scale.factors") |> _[image_scalefactors]
              scale_imgs <- scale_imgs / scale_fct

              # eg image part
              slot(seu_obj, "images")[[fov_section]] |>
                slot("image") |> str()
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
            # sanity, if number of genes doesn't overlap, subset SFE to genes present in mols
            genes_mols <- unique(mols$ID)
            if (!all(genes_mols %in% rownames(sfe))) {
              gene_indx <-
                which(rownames(sfe) %in% genes_mols)
              # genes/features that are removed
              genes_rm <- rownames(sfe)[-gene_indx]
              warning(">>> Total number of genes in SFE object does NOT overlap with `transcript` coordinates genes", "\n",
                      ">>> Subsetting object and `transcript` coordinates to keep only matching genes", "\n",
                      "Total of ", length(genes_rm), " genes are removed")
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

          # add dimensionality reduction if present ----
          dimRed_names <- slot(seu_obj, "reductions") |> names()
          if (!is.null(dimRed_names)) {
            message(">>> Adding Dimensionality Reduction: ", paste0(dimRed_names, collapse = ", "))
            dimRed <- slot(seu_obj, "reductions")
            # sanity, make sure cell ids match between dim. reductions and sfe object
            if (any(!rownames(dimRed[[1]]) %in% colnames(sfe)))
              dimRed <-
                sapply(names(dimRed),
                       function(i) base::subset(dimRed[[i]], colnames(sfe)))
            for (i in seq(dimRed_names))
              reducedDim(x = sfe, type = toupper(dimRed_names[i])) <-
                dimRed[[dimRed_names[i]]] |> slot("cell.embeddings")
            # PCA specifics
            if (grep("PCA", reducedDimNames(sfe)))
              # change naming, eg PC_1 to PC1
              for (i in c("cell.embeddings", "feature.loadings"))
                attr(slot(dimRed[["pca"]], i), "dimnames")[[2]] <-
                  gsub("_", replacement = "", dimnames(dimRed[["pca"]])[[2]])
              # add PC variance
              tot.var <- slot(dimRed[["pca"]], "misc") |> unlist()
              varExplained <- slot(dimRed[["pca"]], "stdev")^2
              percentVar <- (varExplained / tot.var) * 100
              attr(reducedDim(sfe, type = "PCA"), "varExplained") <- varExplained
              attr(reducedDim(sfe, type = "PCA"), "percentVar") <- percentVar
              # add PC loadings
              attr(reducedDim(sfe, type = "PCA"), "rotation") <-
              dimRed[["pca"]] |> slot("feature.loadings")
          }

          # TODO (alternatively): use `MultiAssayExperiment` instead? ----
          # Currently using `altExp` if > 1 Seurat Assays are present ----
          if (length(assays_name) > 1) {
            # keep assays minus default one
            assays_name <- setdiff(assays_name, .DefaultAssay(seu_obj))
            # sanity, check if assays cell ids match those in fov_seu
            cells_passed <-
              lapply(assays_name, function(i)
                Seurat::GetAssay(seu_obj, i) |>
                  SeuratObject::Cells() |>
                  match(x = _, cell_ids_fov) |>
                  na.omit() |> any()) |> unlist()
            if (all(cells_passed)) {
              message(">>> Adding Seurat Assay(s) as Alternative Experiment(s): ",
                      paste0("\n", assays_name))
              # loop and add multiple assays
              for (a in seq(assays_name)) {
                # get slot or layer names
                slot_names <- .GetSlotNames(object_seu = seu_obj,
                                            assay_seu = assays_name[a],
                                            fov_number = fov_section)
                # prepare assays
                assays_add <-
                  list(counts = .GetCounts(seu_obj,
                                           cells = cell_ids_fov,
                                           if (!slot_names[["slots_ok"]] && length(slot_names) == 1)
                                             slot_names
                                           else if (!slot_names[["slots_ok"]] && length(slot_names[["slot_n"]]) == 2)
                                             slot_names[["slot_n"]][1]
                                           else "counts",
                                           assay = assays_name[a]),
                       logcounts = .GetCounts(seu_obj,
                                              cells = cell_ids_fov,
                                              if (!slot_names[["slots_ok"]] && length(slot_names[["slot_n"]]) == 2)
                                                slot_names[2]
                                              else "data",
                                              assay = assays_name[a]),
                       scaledata = .GetCounts(seu_obj,
                                              cells = cell_ids_fov,
                                              "scale.data",
                                              assay = assays_name[a]))
                # remove empty elements
                assays_add <- assays_add[lapply(assays_add, length) > 0]
                altExp(sfe, assays_name[a]) <-
                  SingleCellExperiment(assays = assays_add) |>
                  toSpatialExperiment(sample_id = gsub("_|-|[.]", "", x = fov_names[fov_section])) |>
                  toSpatialFeatureExperiment(x = _,
                                             spatialCoordsNames = c("x", "y"),
                                             colGeometries = list(centroids = colGeometry(sfe)),
                                             unit = unit(sfe))
              }
            }
          }
          return(sfe)
        })

      # remove this afterwards
      #obj_list_test <<- obj_list

      if (length(obj_list) > 1) {
        message(">>> Combining ", length(obj_list),
                " SFE object(s) with unique `sample_id`")

        # TODO: issue when one sfe has: ----
          # assays(2): counts logcounts
          # and assays(1): counts
          # Error in .bind_Assays_objects(objects, along.cols = TRUE) :
          #  the objects to bind must have the same number of assays
          #

        sfe <- do.call(cbind, obj_list)
        return(sfe)

        # TODO (alternatively) return list of objects? (when multiple samples/FOVs) ----
        #return(obj_list)
      } else if (length(obj_list) == 1) {
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
