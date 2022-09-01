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
groups <- read.csv(snakemake@params[["sample_sheet"]],
    sep = "\t", header = TRUE)
model_name <- snakemake@params[["groups"]]
sel_model <- groups %>% select(model_name)
model <- sel_model[!(is.na(sel_model) | sel_model == ""), ]
anno <- data.frame(model, row.names = colnames(sleuth_file))
pdf(file = snakemake@output[[1]], height = 10, width = 10)
pheatmap(sleuth_file[selectedgenes, ], annotation = anno, scale = "row")
dev.off()