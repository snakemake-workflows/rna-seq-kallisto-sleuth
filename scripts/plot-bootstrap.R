suppressMessages({
  library("sleuth")
})

so <- sleuth_load(snakemake@input[[1]])

pdf(file = snakemake@output[[1]])
plot_bootstrap(so, snakemake@wildcards[["transcript"]], units = "tpm")
dev.off()
