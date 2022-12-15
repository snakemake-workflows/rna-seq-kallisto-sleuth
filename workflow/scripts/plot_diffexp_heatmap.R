log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(pheatmap)
library(dplyr)
library(tidyr)
#Reading the sleuth log count matrix file
sleuth_file <- read.csv(snakemake@input[["logcountmatrix_file"]],
    sep = "\t", header = TRUE)
#Adding gene name to corresponding transcript id from sleuth file
rownames(sleuth_file) <- paste(sleuth_file$gene,
    ":", sleuth_file$transcript)
sleuth_file$gene <- NULL
sleuth_file$transcript <- NULL

#Getting top variable genes
vargenes <-
    apply(sleuth_file, 1, var)
#Selecting top 50 variable genes
selectedgenes <-
    names(vargenes[order(vargenes, decreasing = TRUE)][1:50])
pdf(file = snakemake@output[[1]], height = 10, width = 10)
pheatmap(sleuth_file[selectedgenes, ], scale = "row")
dev.off()