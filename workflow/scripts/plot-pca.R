source(snakemake@input[["pca_src"]])

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("ggpubr")

#principal components
pc <- 4

# Load data
so <- sleuth_load(snakemake@input[["rds"]])
covariate_column = snakemake@wildcards[["covariate"]]

# Delete NA values
if (snakemake@params[["exclude_nas"]]) {
  so$sample_to_covariates <- subset(so$sample_to_covariates, !is.na(so$sample_to_covariates[[covariate_column]]))
}

plot_pca(so, color_by = covariate_column)
ggsave(snakemake@output[["pca"]], width=14)

# plot pc variance
plot_pc_variance(so, use_filtered = TRUE, units = "est_counts", pca_number = pc)
ggsave(snakemake@output[["pc_var"]], width=14)

# plot loadings
pc_loading_plots <- list()
for(i in 1:pc) {
  pc_loading_plots[[paste0("pc", i, "_loading")]] <- plot_loadings(so, 
                                                                   scale = TRUE, 
                                                                   pc_input = i,
                                                                   pc_count = 10, 
                                                                   units = "est_counts") +
    ggtitle(paste0(covariate_column, ": plot loadings for principal component ", i)) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
}
plots_loading <- ggarrange(plotlist = pc_loading_plots, ncol = 1, nrow = 1, common.legend = TRUE)
ggexport(plots_loading, filename = snakemake@output[["loadings"]], width = 14)