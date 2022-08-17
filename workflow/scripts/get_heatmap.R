log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(pheatmap)
library(dplyr)
library(tidyr)

sleuth_file <- read.csv(snakemake@input[["Sleuth_logcountmatrix_file"]],
    sep = "\t", header = TRUE)
rownames(sleuth_file) <- paste(sleuth_file$transcript,
    "_", sleuth_file$gene)
sleuth_file$gene <- NULL
sleuth_file$transcript <- NULL

vargenes <-
    apply(sleuth_file, 1, var)
selectedgenes <-
    names(vargenes[order(vargenes, decreasing = TRUE)][1:50])
png(snakemake@output[["png"]], w = 908, h = 953, pointsize = 200)
pheatmap(sleuth_file[selectedgenes, ], scale = "row")
dev.off()