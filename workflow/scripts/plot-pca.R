log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")

so <- sleuth_load(snakemake@input[[1]])
pdf(file = snakemake@output[[1]])
plot_pca(so, color_by = snakemake@wildcards[["covariate"]], units = "tpm", text_labels = TRUE)
dev.off()
