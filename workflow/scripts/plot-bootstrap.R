log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tidyverse")
library("sleuth")

so <- sleuth_load(snakemake@input[["so"]])

top_n <- -strtoi(snakemake@params["top_n"])

results <- read_tsv(snakemake@input[["transcripts"]])

top_genes <- results %>%
   filter(qval <= snakemake@params[["fdr"]]) %>%
   top_n(top_n, qval) %>%
   dplyr::select(ext_gene) %>%
   drop_na() %>%
   distinct(ext_gene)


if (snakemake@params[["genes"]]$activate == TRUE) {
    gene_table <- read.csv(snakemake@params[["genes"]]$genelist,
        sep = "\t")
    names(gene_table) <- c("ext_gene")
    genes_of_interest <- tibble(gene_table) %>%
        distinct(ext_gene)
} else {
    # "genes" is null, if the list provided in config.yaml is empty
    genes_of_interest <- tibble(ext_gene = character())
}
genes <- full_join(top_genes, genes_of_interest, by = "ext_gene") %>%
    add_row(ext_gene = "Custom") %>%
    pull(ext_gene)

dir.create(snakemake@output[[1]])

for (gene in genes) {
    transcripts <- results %>%
        filter(ext_gene == gene) %>%
        drop_na() %>%
        pull(target_id)

    if (length(transcripts > 0)) {
        for (transcript in transcripts) {
            plot_bootstrap(so, transcript, color_by = snakemake@params[["color_by"]], units = "est_counts")
            ggsave(file = str_c(snakemake@output[[1]], "/", gene, ".", transcript, ".", snakemake@wildcards[["model"]], ".bootstrap.pdf"))
        }
    }
}
