log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("ggpubr")

so <- sleuth_load(snakemake@input[[1]])
pc12 <- plot_pca(so, color_by = snakemake@wildcards[["covariate"]], units = "est_counts", text_labels = TRUE)
pc34 <- plot_pca(so, pc_x = 3, pc_y = 4, color_by = snakemake@wildcards[["covariate"]], units = "est_counts", text_labels = TRUE)

ggarrange(pc12, pc34, ncol = 2, common.legend = TRUE)
ggsave(snakemake@output[[1]], width=14)
