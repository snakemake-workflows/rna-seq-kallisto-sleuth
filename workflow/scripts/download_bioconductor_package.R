log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(BiocManager)

BiocManager::install("org.Hs.eg.db", lib = snakemake@params[["lib"]])
