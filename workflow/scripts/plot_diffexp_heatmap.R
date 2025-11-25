log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(tidyverse)
library(ggalign)

# Reading the sleuth log count matrix file
log_counts <- read_tsv(
    snakemake@input[["logcountmatrix_file"]]
)

# Check the mode of the rule's execution, topn vs. predefined gene list
if (snakemake@wildcards[["gene_list"]] == "topn") {
    selected_genes <- log_counts |>
        pivot_longer(
            !c(transcript, gene),
            names_to = "sample",
            values_to = "log_count"
        ) |>
        summarize(
            .by = c(gene, transcript),
            variance = var(log_count)
        ) |>
        left_join(
            log_counts,
            by = join_by(gene, transcript)
        ) |>
        slice_max(variance, n = 50) |>
        select(-variance)
} else {
    # Adding gene list to the variable
    predef_gene_list <-
        read_tsv(
            snakemake@input[["predef_gene_list"]],
            col_names = c("gene")
        )
    selected_genes <-
        log_counts |>
        inner_join(
            predef_gene_list,
            by = join_by(gene)
        )
}

matrix <- selected_genes |>
    mutate(transcript = str_pad(transcript, 18, "right")) |>
    replace_na(list(gene = "   ")) |>
    mutate(t = str_c(gene, " - ", transcript)) |>
    select(-gene, -transcript) |>
    column_to_rownames(var = "t")

heatmap <- ggheatmap(matrix) +
    geom_tile(aes(fill = value), color = "grey") +
    theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
    anno_top(size = unit(3, "cm")) +
    align_dendro() +
    theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
    ) +
    anno_right(size = unit(3, "cm")) +
    align_dendro() +
    theme(
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
    )

ggsave(
    snakemake@output[["diffexp_heatmap"]],
    heatmap,
    width = length(colnames(matrix)) * 0.5 + 6,
    height = length(rownames(matrix)) * 0.2 + 3,
    limitsize = FALSE
)
