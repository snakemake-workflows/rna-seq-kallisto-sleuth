log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")

so <- sleuth_load(snakemake@input[[1]])

pdf(file = snakemake@output[[1]])
plot_fld(so, paste0(snakemake@wildcards[["sample"]], "-", snakemake@wildcards[["unit"]]))
dev.off()
