log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(pheatmap)
library(dplyr)
library(tidyr)

sleuth_file <- read.csv(snakemake@input[["logcountmatrix_file"]],
    sep = ",", header = TRUE)
rownames(sleuth_file) <- paste(kallsito_file$transcript,
    "_", sleuth_file$gene)
sleuth_file$gene < -NULL
sleuth_file$transcript <- NULL

predefine_genelist <- read.table(snakemake@input[["predef_genelist"]],
   sep = " ", header = TRUE)

selectedgenes <-
    sleuth_file %>% filter(rownames(sleuth_file)
     %in% predefine_genelist$Selected_genes)

print(selectedgenes)
png(snakemake@output[["predefgene_png"]], w = 1500, h = 1500, pointsize = 200)
pheatmap(selectedgenes, scale = "row")
dev.off()