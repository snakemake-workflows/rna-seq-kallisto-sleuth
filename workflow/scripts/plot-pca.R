log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("ggpubr")

so <- sleuth_load(snakemake@input[[1]])

# number of principal components for pca plot, variance-plot and loadings plot
n_pc <- snakemake@params[["n_pc"]]
n_pca <- n_pc

if(n_pca%%2 != 0){
  n_pca <- n_pca + 1
}

# plot pca
pca_plot_list <- list()
i <- 1
while(i < n_pca){
  print(paste0("pca_pc", i, i + 1))
  pca_plot_list[[paste0("pca_pc", i, i + 1)]] <-
    plot_pca(so, pc_x = i, pc_y = i + 1, color_by = snakemake@wildcards[["covariate"]],
             use_filtered = TRUE, units = "est_counts", text_labels = TRUE)+
    ggtitle(paste0(snakemake@wildcards[["covariate"]], ": pca plot for principal components ", i, " vs. ", i+1))+
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  i <- i + 2
}

pca_plots <- ggarrange(plotlist = pca_plot_list, ncol = 1, nrow = 1, common.legend = TRUE)
ggexport(pca_plots, filename = snakemake@output[["pca"]])

# plot pc variance
plot_pc_variance(so, use_filtered = TRUE, units = "est_counts", pca_number = n_pc)
ggsave(snakemake@output[["pc_var"]], width=14)

# plot loadings
pc_loading_plots <- list()
for(i in 1:n_pc) {
  pc_loading_plots[[paste0("pc", i, "_loading")]] <- plot_loadings(so, scale = TRUE, pc_input = i,
                                                                   pc_count = n_pc, units = "est_counts")+
    ggtitle(paste0(snakemake@wildcards[["covariate"]], ": plot loadings for principal component ", i))+
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
}
plots_loading <- ggarrange(plotlist = pc_loading_plots, ncol = 2, nrow = 1, common.legend = TRUE)
ggexport(plots_loading, filename = snakemake@output[["loadings"]], width = 14)