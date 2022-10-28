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
#Adding corresponding sample groups
groups <- read.csv(snakemake@params[["sample_sheet"]],
    sep = "\t", header = TRUE)
model_name <- snakemake@params[["groups"]]
sel_model <- groups %>% select(model_name)
model <- sel_model[!(is.na(sel_model) | sel_model == ""), ]
anno <- data.frame(model, row.names = colnames(sleuth_file))
colnames(anno) <- snakemake@params[["groups"]]
#Plotting the heatmap
pdf(file = snakemake@output[[1]], height = 10, width = 10)

pheatmap(sleuth_file[selectedgenes, ], annotation = anno, scale = "row")

dev.off()