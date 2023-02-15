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
if (snakemake@params[["predef_genelist"]]$activate == TRUE) {
    predefine_genelist <-
        read.table(snakemake@params[["predef_genelist"]]$genelist,
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
    if (all(selectedgenes == 0)) {
        txt <- "cannot plot, all values are zero"
        pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
        plot.new()
        text(.5, .5, txt, font = 2, cex = 1.5)
        dev.off()
    # If number of transcripts are more than 100,
    # Since for RNA-seq multiple transcripts may present for single gene.
    } else if (nrow(selectedgenes) > 100) {
        pdf_height <- round(nrow(selectedgenes) / 7)
        pdf(snakemake@output[["diffexp_heatmap"]],
            height = pdf_height, width = ncol(selectedgenes)
        )
        pheatmap(selectedgenes, selectedgenes,
            cluster_rows = FALSE,  display_numbers = TRUE,
                cellheight = 10, scale = "row"
        )
        dev.off()
    } else {
        fontsize_row <- 10 - nrow(selectedgenes) / 15
        pdf(snakemake@output[["diffexp_heatmap"]])
        pheatmap(selectedgenes, selectedgenes,
            cluster_rows = FALSE, fontsize_row = fontsize_row, scale = "row"
        )
        dev.off()
    }
} else {
    # Adding gene name to corresponding transcript id from sleuth file
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
    selectedgenes <-
        names(vargenes[order(vargenes, decreasing = TRUE)][1:50])
    pdf(file = snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
    pheatmap(sleuth_file[selectedgenes, ], scale = "row")
    dev.off()
}
