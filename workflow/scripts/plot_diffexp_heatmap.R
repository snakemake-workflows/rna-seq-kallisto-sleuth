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
# Check the config file if it contains pre-defined gene list

if (snakemake@wildcards[["mode"]] == "topn") {
    # Replacing the rownames with transcript Id and gene name
    rownames(sleuth_file) <- paste(
        sleuth_file$transcript,
        ":", sleuth_file$gene
    )
    sleuth_file$gene <- NULL
    sleuth_file$transcript <- NULL
    # Getting top variable genes
    vargenes <-
        apply(sleuth_file, 1, var)
    # Selecting top 50 variable genes
    selectedgenes <- sleuth_file[row.names(sleuth_file)
    %in% names(vargenes[order(vargenes, decreasing = TRUE)][1:50]), ]
    # If the config file contains pre-defined gene list
} else if (snakemake@wildcards[["mode"]] == "predefined") {
    # Adding gene list to the variable
    predefine_genelist <-
        read.table(snakemake@input[["predef_genelist"]],
            sep = "\t"
        )
    selectedgenes <-
        sleuth_file %>% filter(sleuth_file$gene
            %in% predefine_genelist$V1)
    rownames(selectedgenes) <- paste(
        selectedgenes$transcript,
        ":", selectedgenes$gene
    )
    selectedgenes$gene <- NULL
    selectedgenes$transcript <- NULL
}
# Checks if selectedgenes variable contain genes that have expression values
if (all(selectedgenes == 0)) {
    txt <- "cannot plot, all values are zero"
    pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
    plot.new()
    text(.5, .5, txt, font = 2, cex = 1.5)
    dev.off()
    # If number of transcripts is more than 100,
    # Since for RNA-seq multiple transcripts may present for a single gene.
} else if (nrow(selectedgenes) > 100) {
    pdf_height <- round(nrow(selectedgenes) / 5)
    pdf(snakemake@output[["diffexp_heatmap"]],
        height = pdf_height, width = ncol(selectedgenes) * 2
    )
    pheatmap(selectedgenes, selectedgenes,
        cluster_rows = FALSE, display_numbers = TRUE,
        cellheight = 10, scale = "row"
    )
    dev.off()
} else {
    pdf_height <- round(nrow(selectedgenes) / 5)
    pdf(
        file = snakemake@output[["diffexp_heatmap"]],
        height = pdf_height, width = ncol(selectedgenes) * 2
    )
    pheatmap(selectedgenes, scale = "row")
    dev.off()
}