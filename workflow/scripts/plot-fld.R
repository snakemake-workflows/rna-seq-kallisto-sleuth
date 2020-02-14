log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("ggplot2")

so <- sleuth_load(snakemake@input[[1]])

p <- plot_fld(so, snakemake@wildcards[["sample"]])

ggsave(filename = snakemake@output[[1]], plot = p)
