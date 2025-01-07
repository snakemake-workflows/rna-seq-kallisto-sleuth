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
diffexp_file <- read.csv(snakemake@input[["diffexp_file"]],
    sep = "\t", header = TRUE
)

sample_file <- read.csv(snakemake@input[["sample_file"]],
    sep = "\t", header = TRUE
)


# Replacing the rownames with transcript Id and gene name
rownames(sleuth_file) <- paste(
    sleuth_file$transcript,
    ":", sleuth_file$gene
)
sleuth_file$gene <- NULL
sleuth_file$transcript <- NULL

# Compute name of genes with biggest absolute of signed_pi_values
number_genes <- snakemake@params[["number_genes"]]
top_gene_strings <- paste(diffexp_file$target_id, ":", diffexp_file$ext_gene)
selectedgenes <- sleuth_file[row.names(sleuth_file) %in% top_gene_strings[1:number_genes], ]
selectedgenes <- selectedgenes[match(top_gene_strings[1:number_genes], row.names(selectedgenes)), ]



if (all(selectedgenes == 0)) {
    txt <- "Cannot plot, all values are zero"
    pdf(snakemake@output[["diffexp_heatmap"]], height = 10, width = 10)
    plot.new()
    text(.5, .5, txt, font = 2, cex = 1.5)
    dev.off()
} else {
    pdf_height <- max(10, round(nrow(selectedgenes) / 5))  # Minimum-Höhe von 10
    pdf(
        file = snakemake@output[["diffexp_heatmap"]],
        height = pdf_height, width = max(10, ncol(selectedgenes) * 2)
    )

    # Compute group for each sample 
    model <- snakemake@params[["model"]]$primary_variable
    group1 <- strsplit(model, split = '_vs_')[[1]][1]
    group2 <- strsplit(model, split = '_vs_')[[1]][2]
    sample_group <- sample_file %>%
        mutate(
            !!model := ifelse(
            .[[model]] == "+", group1, 
            ifelse(.[[model]] == "-", group2, NA)
            )) %>%
        mutate(sample = paste0("X", sample, ".0.1")) %>% 
        filter(!is.na(.[[model]])) %>%
        select(sample, !!model)
    
    rownames(sample_group) <- sample_group$sample  
    sample_group$sample <- NULL
    pheatmap(selectedgenes, scale = "row", cluster_rows = FALSE, cluster_cols = TRUE, annotation_col = sample_group)
    dev.off()
}