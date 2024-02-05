library(vroom)
library(BiocParallel)
library(tidyverse)
# Experiment-------------
df <- vroom("~/detected_transcripts_small.csv",
            col_select = c(1, global_x:global_z, gene, transcript_id),
            col_types = vroom::cols(...1 = "c"))
df <- vroom("~/detected_transcripts.csv",
            col_select = c(1, global_x:global_z, gene, transcript_id),
            col_types = vroom::cols(...1 = "c"))
genes <- df |>
    select(gene, transcript_id) |>
    distinct()
anyDuplicated(genes$gene)
anyDuplicated(genes$transcript_id)
genes[343,]
genes |> filter(str_detect(gene, "Blank"))
# OK, one transcript for one gene, except for the blanks.
# Shall I care about multiple isoforms? Maybe not yet. I don't see it coming anytime soon.
unique(df$global_z)

df <- df[df$global_z == 3, c("global_x", "global_y", "gene")]
names(df) <- c("x", "y", "group")
# Not as slow as I expected
df_sf <- df2sf(df, geometryType = "MULTIPOINT")
object.size(df_sf) |> format(units = "MB")
# Not as bad as I thought
saveRDS(df_sf, "merfish_liver_spots.rds")

# Make test dataset from FOV1-----------
df <- vroom("~/detected_transcripts.csv",
            col_types = vroom::cols(...1 = "c"))
df <- df |> filter(fov == 0)
# Are different z planes different?
df |>
    filter(gene == "Ldha") |>
    ggplot(aes(global_x, global_y, color = global_z)) +
    geom_point(shape = 3)
# Answer: yes.
df <- column_to_rownames(df, "...1")
nmols <- nrow(df)
df_sample <- df[sample(seq_len(nmols), round(nmols/50)),]
write.csv(df_sample, "inst/extdata/vizgen/detected_transcripts.csv",
          quote = FALSE, row.names = TRUE)

ggplot(df_sample, aes(global_x, global_y, color = global_z)) +
    geom_point(shape = 3)

# Xenium---------
library(arrow)

transcripts <- read_parquet("xenium_skin/transcripts.parquet")
transcripts |>
    filter(fov_name == "G1", between(x_location, 500, 520),
           between(y_location, 3700, 3720)) |>
    ggplot(aes(x_location, y_location, color = z_location)) +
    geom_point(shape = 3) +
    coord_equal()

ggplot(transcripts, aes(qv, fill = cell_id == "UNASSIGNED")) +
    geom_histogram(bins = 50)

# CosMX---------
cosmx_tx <- vroom("cosmx_brain/Quarter Brain/Run5642_S3_Quarter_tx_file.csv")
