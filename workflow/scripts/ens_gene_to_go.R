log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(snakemake@params[["bioc_species_pkg"]], character.only = TRUE)

# provides `tidyverse` and load_bioconductor_package()
source(snakemake@params[["common_src"]])

ens_gene_to_go <-
    AnnotationDbi::select(get(snakemake@params[["bioc_species_pkg"]]),
        keys = keys(get(snakemake@params[["bioc_species_pkg"]]),
            keytype = "ENSEMBL"
        ), columns = c("GO"),
        keytype = "ENSEMBL"
    ) %>%
    # rows with empty ENSEMBL or GO IDs are useless and will lead to trouble below
    drop_na(ENSEMBL, GO) %>%
    dplyr::rename(ensembl_gene_id = ENSEMBL) %>%
    group_by(ensembl_gene_id) %>%
    summarise(go_ids = str_c(GO, collapse = ";"))


write_tsv(ens_gene_to_go, path = snakemake@output[[1]], col_names = FALSE)
