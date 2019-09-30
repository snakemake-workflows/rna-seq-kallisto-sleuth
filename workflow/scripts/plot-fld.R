suppressMessages({
  library("sleuth")
})

so <- sleuth_load(snakemake@input[[1]])

pdf(file = snakemake@output[[1]])
plot_fld(so, paste0(snakemake@wildcards[["sample"]], "-", snakemake@wildcards[["unit"]]))
dev.off()
