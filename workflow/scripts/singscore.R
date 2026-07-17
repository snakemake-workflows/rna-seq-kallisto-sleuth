log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(singscore)
library(GSEABase)
library(ComplexHeatmap)


#### get data for calculation 
xpr_data  <- read_tsv(snakemake@input[["gene_counts"]])
smpl_data <- read_tsv(snakemake@input[["samples"]])
gene_set  <- getGmt(snakemake@input[["gene_sets"]])
transcript_index <- read_tsv(snakemake@input[["transcript_ref"]])

group <- snakemake@wildcards[["gene_set_group"]]

#### definition of aes for the plot 
color_aes <- snakemake@params[["color_aes"]]
shape_aes <- snakemake@params[["shape_aes"]]

#### calculation of singscore 
#### get only mane transcripts 
trans_mane  <- transcript_index %>%
  filter(canonical == "TRUE",
         !is.na(ens_gene))

xpr_matrix <- xpr_data %>%
  mutate(transcript = sub("\\..*", "", transcript),
         row_flt = rowSums(across(where(is.numeric)))) %>%
  arrange(desc(row_flt)) %>%
  filter(transcript %in% trans_mane$target_id) %>%
  mutate(ens_gene = trans_mane$ens_gene[match(transcript, trans_mane$target_id)]) %>%
  distinct(ens_gene, .keep_all = TRUE) %>%
  column_to_rownames(var = "ens_gene") %>%
  dplyr::select(-c("transcript", "row_flt", "gene")) %>%
  as.matrix()

data_ranked      <- rankGenes(xpr_matrix)


#### get genesets 
all_sets <- names(gene_set)

up_set   <- all_sets[grepl("up", all_sets, ignore.case = TRUE)]
down_set <- all_sets[grepl("down", all_sets, ignore.case = TRUE)]

singscore_output <- if (!is.null(down_set)) {
  simpleScore(data_ranked, upSet = gene_set[[up_set]], downSet = gene_set[[down_set]])
} else {
  simpleScore(data_ranked, upSet = gene_set[[up_set]])
}

score_df <- data.frame(
  sample = colnames(xpr_matrix),
  score  = singscore_output$TotalScore) %>%
  mutate(sample = sub("t-.*", "t", sample)) %>%
  left_join(smpl_data, by = "sample") %>%
  distinct(sample, .keep_all = TRUE) %>%
  arrange(score)

# ---- plot, sized dynamically to the number of samples ----
n_samples   <- nrow(score_df)
plot_width  <- max(6, n_samples * 0.15)
plot_height <- 4

tiff(
  filename    = snakemake@output[[1]],
  res         = 150,
  compression = "lzw",
  units       = "in",
  width       = plot_width,
  height      = plot_height,
  pointsize   = 16
)

ggplot(score_df, aes(x = sample, y = score,
                     color = .data[[color_aes]],
                     shape = .data[[shape_aes]])) +
  geom_point(size = 12 / .pt) +
  scale_x_discrete(limits = score_df$sample) +
  labs(
    title = paste(group, "Activity"),
    y     = "SingScore"
  ) +
  theme(
    panel.background     = element_blank(),
    axis.line             = element_line(color = "#000000", linewidth = 1),
    axis.title.x          = element_blank(),
    axis.text.x           = element_text(angle = 270, hjust = .5, vjust = .5,
                                         size = 36 / .pt, color = "#000000"),
    axis.title.y          = element_text(hjust = .5, vjust = .5, size = 52 / .pt),
    axis.text.y            = element_text(size = 36 / .pt, color = "#000000"),
    legend.title           = element_blank(),
    legend.text            = element_text(color = "#000000", size = 36 / .pt),
    legend.justification   = "top"
  )

dev.off()

sink(type = "message")
sink(type = "output")