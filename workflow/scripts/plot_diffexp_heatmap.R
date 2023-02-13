log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(pheatmap)
library(dplyr)
library(tidyr)
# Reading the sleuth log count matrix file
sleuth_file <- read.csv(snakemake@input[["logcountmatrix_file"]],
    sep = "\t", header = TRUE
)
print(length(snakemake@params[["predef_genelist"]]))
print(snakemake@params[["predef_genelist"]]$activate)
# Getting pre-defined genelist file
if (snakemake@params[["predef_genelist"]]$activate == TRUE &&
    length(snakemake@params[["is_3prime_experiment"]]) == "true") {
    predefine_genelist <- read.table(snakemake@params[["predef_genelist"]],
        sep = "\t"
    )
    selectedgenes <-
        sleuth_file %>% filter(sleuth_file$gene
            %in% predefine_genelist$V1)
    rownames(selectedgenes) <- paste(
        selectedgenes$gene,
        ":", selectedgenes$transcript
    )
    selectedgenes$gene <- NULL
    selectedgenes$transcript <- NULL
    if (all(selectedgenes == 0)) {
        txt <- "cannot plot, all values are zero"
        pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
        plot.new()
        text(.5, .5, txt, font = 2, cex = 1.5)
        dev.off()
    } else {
        pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
        pheatmap(selectedgenes, scale = "row")
        dev.off()
    }
} else if (length(snakemake@params[["is_3prime_experiment"]]) == "true") {
    # Adding gene name to corresponding transcript id from sleuth file
    rownames(sleuth_file) <- paste(
        sleuth_file$gene,
        ":", sleuth_file$transcript
    )
    sleuth_file$gene <- NULL
    sleuth_file$transcript <- NULL
    # Getting top variable genes
    vargenes <-
        apply(sleuth_file, 1, var)
    # Selecting top 50 variable genes
    selectedgenes <-
        names(vargenes[order(vargenes, decreasing = TRUE)][1:50])
    pdf(file = snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
    pheatmap(sleuth_file[selectedgenes, ], scale = "row")
    dev.off()
} else if (snakemake@params[["predef_genelist"]]$activate == TRUE) {
    predefine_genelist <- read.table(snakemake@params[["predef_genelist"]],
        sep = "\t"
    )
    selectedgenes <-
        sleuth_file %>% filter(sleuth_file$gene
            %in% predefine_genelist$V1)
    rownames(selectedgenes) <- paste(
        selectedgenes$transcript
    )
    selectedgenes$gene <- NULL
    selectedgenes$transcript <- NULL
    print(selectedgenes)
    if (all(selectedgenes == 0)) {
        txt <- "cannot plot, all values are zero"
        pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
        plot.new()
        text(.5, .5, txt, font = 2, cex = 1.5)
        dev.off()
    } else {
        pdf(snakemake@output[["diffexp_heatmap"]])
        pheatmap(selectedgenes, cluster_rows = FALSE, scale = "row")
        dev.off()
    }
} else {
    rownames(sleuth_file) <- paste(
        sleuth_file$transcript
    )
    sleuth_file$gene <- NULL
    sleuth_file$transcript <- NULL

    vargenes <-
        apply(sleuth_file, 1, var)
    selectedgenes <-
        names(vargenes[order(vargenes, decreasing = TRUE)][1:50])
    if (all(selectedgenes == 0)) {
        txt <- "cannot plot, all values are zero"
        pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
        plot.new()
        text(.5, .5, txt, font = 2, cex = 1.5)
        dev.off()
    } else {
        pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
        pheatmap(sleuth_file[selectedgenes, ], scale = "row")
        dev.off()
    }
}
