log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("ggpubr")
library("sleuth")
library("cowplot")
rlang::global_entrace()

# Transform to df (Taken from sleuth https://github.com/pachterlab/sleuth/blob/20e493fb98cb1b65a508e6634db430d5d3848567/R/sleuth.R#L928)
spread_abundance_by <- function(abund, var, which_order) {
  abund <- data.table::as.data.table(abund)
  var_spread <- data.table::dcast(abund, target_id ~ sample, value.var = var)
  var_spread <- var_spread[order(var_spread$target_id), ]
  var_spread <- as.data.frame(var_spread, stringsAsFactors = FALSE)
  rownames(var_spread) <- var_spread$target_id
  var_spread["target_id"] <- NULL
  result <- as.matrix(var_spread)
  result[, which_order, drop = FALSE]
}

prepare_pca_df <- function(obj, color_by) {
  # Extract data
  mat <- t(spread_abundance_by(obj$obs_norm_filt, "est_counts", obj$sample_to_covariates$sample))
  # Remove zero variance columns
  non_zero_columns <- apply(mat, 2, function(col) any(col != 0))
  mat_cleaned <- mat[, non_zero_columns]

  # Do pca
  pca_res <- prcomp(mat_cleaned, scale = TRUE)

  # Add color information  
  pca_df <- as.data.frame(pca_res$x)
  pca_df$sample <- rownames(pca_df)
  pca_df <- merge(pca_df, obj$sample_to_covariates[, c("sample", color_by)], by = "sample")
  return(pca_df)
}

#principal components
pc <- 4

# Load data
so <- sleuth_load(snakemake@input[["rds"]])
covariate_column = snakemake@wildcards[["covariate"]]

# Delete NA values
if (snakemake@params[["exclude_nas"]]) {
  so$sample_to_covariates <- subset(so$sample_to_covariates, !is.na(so$sample_to_covariates[[covariate_column]]))
}

pca_df <- prepare_pca_df(so, color_by = covariate_column)
write.table(pca_df, file=snakemake@output[["pca"]], quote=FALSE, sep='\t')

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