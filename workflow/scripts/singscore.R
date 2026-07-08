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
gene_sets <- getGmt(snakemake@input[["gene_sets"]])

#### definition of gene set groups for the analysis 
group <- snakemake@wildcards[["gene_set_group"]]
group_cfg <- snakemake@config[["singscore"]][["gene_sets"]][[group]]

up_name   <- group_cfg[["upregulated_gene_set"]]
down_name <- group_cfg[["downregulated_gene_set"]]

up_set   <- gene_sets[[up_name]]
down_set <- if (!is.null(down_name)) gene_sets[[down_name]] else NULL


#### definition of aes for the plot 
color_aes <- snakemake@params[["color_aes"]]
shape_aes <- snakemake@params[["shape_aes"]]

#### calculation of singscore 
xpr_matrix <- xpr_data %>%
  column_to_rownames(var = "gene_ext") %>%
  as.matrix()
data_ranked      <- rankGenes(xpr_matrix)

singscore_output <- if (!is.null(down_set)) {
  simpleScore(data_ranked, upSet = up_set, downSet = down_set)
} else {
  simpleScore(data_ranked, upSet = up_set)
}

score_df <- data.frame(
  sample = colnames(xpr_matrix),
  score  = singscore_output$TotalScore) %>%
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


