log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")
library("ggpubr")

so <- sleuth_load(snakemake@input[[1]])
plot_list <- list()

# combinations between samples without repetition 

# sample_matrix <- t(combn(so$sample_to_covariates$sample, 2))
# combinations <- as.numeric(row(sample_matrix))
# 
# for(i in combinations) {
#   sample1 <- sample_matrix[i, 1]
#   sample2 <- sample_matrix[i, 2]
#   
#   plot_list[[i]] <- plot_scatter(so, 
#                                  sample_x = sample1, 
#                                  sample_y = sample2, 
#                                  use_filtered = TRUE, 
#                                  units = "est_counts", 
#                                  point_alpha = 0.2, 
#                                  xy_line = TRUE, 
#                                  xy_line_color = "red") +
#     labs(x = sample1, y = sample2)
#   }
  

# all combinations between samples (with repitition)

samples <- so$sample_to_covariates$sample
sample_matrix <- crossing(x = samples, y = samples) # creates all combinations between samples
combinations <- as.numeric(row(sample_matrix))

for(i in combinations) {
  sample1 <- sample_matrix[i, 1]
  sample2 <- sample_matrix[i, 2]
  if(sample1 != sample2) {
    plot_list[[i]] <- plot_scatter(so,
                                   sample_x = sample1,
                                   sample_y = sample2,
                                   use_filtered = TRUE,
                                   units = "est_counts",
                                   point_alpha = 0.2,
                                   xy_line = TRUE,
                                   xy_line_color = "red") +
      labs(x = sample1, y = sample2)

  } else {
    x_val <- y_val <- c(0,0)
    test <- tibble(x_val, y_val)
    plot_list[[i]] <- ggplot(test, aes(x = x_val, y = y_val)) +
      theme_void() +
      geom_text(label = "") +
      annotate("text", label = sample1, x=0, y=0, colour = "blue", fontface = "bold")
  }
}

all_plots <- ggarrange(plotlist = plot_list)
plot_figure <- annotate_figure(all_plots, 
                               top = text_grob("log transformed scatter plots of different samples",
                                               color = "black", 
                                               face = "bold", 
                                               size = 14))
ggsave(snakemake@output[[1]], plot_figure)
