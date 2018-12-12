suppressMessages({
  library("sleuth")
})

so <- sleuth_load(snakemake@input[[1]])

svg(file = snakemake@output[[1]])
plot_pca(so, color_by = snakemake@wildcards[["covariate"]])
dev.off()
