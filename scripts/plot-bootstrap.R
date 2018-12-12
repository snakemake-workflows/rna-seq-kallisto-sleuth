suppressMessages({
  library("sleuth")
})

so <- sleuth_load(snakemake@input[[1]])

svg(file = snakemake@output[[1]])
plot_bootstrap(so, snakemake@wildcards[["transcript"]], units = "est_counts", color_by = "condition")
dev.off()
