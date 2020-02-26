log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")
library("ggpubr")

so <- sleuth_load(snakemake@input[[1]])

sample_matrix <- t(combn(so$sample_to_covariates$sample, 2))
combinations <- as.numeric(row(sample_matrix))

plot_list <- list()

for(i in combinations) {
  sample1 <- sample_matrix[i, 1]
  sample2 <- sample_matrix[i, 2]
  plot_list[[i]] <- plot_scatter(so, sample_x = sample1, sample_y = sample2, use_filtered = TRUE,
               units = "est_counts", offset = 1, point_alpha = 0.2, xy_line = TRUE,
               xy_line_color = "red", trans = "log", xlim = NULL, ylim = NULL)+ 
    labs(x = sample1, y = sample2)
}

all_plots <- ggarrange(plotlist = plot_list)
plot_figure <- annotate_figure(all_plots, 
                top = text_grob("log transformed scatter plots of different samples",
                                color = "black", face = "bold", size = 14)
                )
ggsave(snakemake@output[[1]], plot_figure)