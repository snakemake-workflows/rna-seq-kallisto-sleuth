log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tidyverse")
library("sleuth")

so <- sleuth_load(snakemake@input[["so"]])

top_n <- -strtoi(snakemake@params["top_n"])

transcripts <- read_tsv(snakemake@input[["transcripts"]])

top_genes <- transcripts |>
    filter(qval <= snakemake@params[["fdr"]]) |>
    top_n(top_n, qval) |>
    dplyr::select(ext_gene) |>
    drop_na() |>
    distinct(ext_gene)

if (snakemake@params[["genes_of_interest"]][["activate"]] == TRUE) {
    genes_of_interest <- read_tsv(
        unlist(unname(snakemake@params[["genes_of_interest"]][["gene_lists"]])),
        col_names = c("ext_gene")
    ) |>
        unique()
} else {
    # Create empty tibble when genes_of_interest is not activated
    genes_of_interest <- tibble(ext_gene = character())
}

genes <- full_join(top_genes, genes_of_interest, by = "ext_gene") |>
    add_row(ext_gene = "Custom")

dir.create(snakemake@output[[1]])

# sleuth seems to do some internal filtering, and bootstrap values are only
# available for filtered transcripts
available_transcripts <- so$obs_norm_filt |>
    select(target_id)

transcripts_of_interest <- inner_join(
    transcripts,
    genes,
    by = "ext_gene"
) |>
    inner_join(
        available_transcripts,
        by = "target_id"
    ) |>
    pull(target_id)


if (length(transcripts_of_interest > 0)) {
    for (transcript in transcripts_of_interest) {
        gene = transcripts |> filter(target_id == transcript) |> pull(ext_gene)
        plot_bootstrap(
            so,
            transcript,
            color_by = snakemake@params[["color_by"]],
            units = "est_counts"
        )
        ggsave(
            file = str_c(
                snakemake@output[[1]],
                "/",
                gene,
                ".",
                transcript,
                ".",
                snakemake@wildcards[["model"]],
                ".bootstrap.pdf"
            )
        )
    }
}
