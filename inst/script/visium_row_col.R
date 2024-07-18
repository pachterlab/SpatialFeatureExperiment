# I installed Space Ranger in my home directory
spacerange_dir <- "~/tmp/spaceranger-3.0.1/lib/python/cellranger/barcodes"

get_data <- function(file_path) {
    data <- read.delim(file_path, header = FALSE, col.names = c("barcode", "col", "row"))
    return(data)
}

visium_row_col <- get_data(paste0(spacerange_dir, "/visium-v1_coordinates.txt"))
visium_row_col_v1 <- get_data(paste0(spacerange_dir,"/visium-v1_coordinates.txt"))
visium_row_col_v2 <- get_data(paste0(spacerange_dir,"/visium-v2_coordinates.txt"))
visium_row_col_v3 <- get_data(paste0(spacerange_dir,"/visium-v3_coordinates.txt"))
visium_row_col_v4 <- get_data(paste0(spacerange_dir,"/visium-v4_coordinates.txt"))
visium_row_col_v5 <- get_data(paste0(spacerange_dir,"/visium-v5_coordinates.txt"))

usethis::use_data(visium_row_col, overwrite = TRUE)
usethis::use_data(visium_row_col_v1, overwrite = TRUE)
usethis::use_data(visium_row_col_v2, overwrite = TRUE)
usethis::use_data(visium_row_col_v3, overwrite = TRUE)
usethis::use_data(visium_row_col_v4, overwrite = TRUE)
usethis::use_data(visium_row_col_v5, overwrite = TRUE)

