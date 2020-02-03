log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# provides `tidyverse` and load_bioconductor_package()
source( file.path(snakemake@scriptdir, 'common.R') )

pkg <- snakemake@params[["bioc_pkg"]]
load_bioconductor_package(snakemake@input[["species_anno"]], pkg)

ens_gene_to_go <- AnnotationDbi::select(  get(pkg),
                                          keys=keys(get(pkg), keytype="ENSEMBL"),
                                          columns=c("GO"),
                                          keytype="ENSEMBL"
                                          ) %>%
                  # rows with empty ENSEMBL or GO IDs are useless and will lead to trouble below
                  drop_na(ENSEMBL,GO) %>%
                  dplyr::rename(ensembl_gene_id = ENSEMBL) %>%
                  group_by(ensembl_gene_id) %>%
                  summarise(go_ids = str_c(GO, collapse = ";"))


write_tsv(ens_gene_to_go, path = snakemake@output[[1]], col_names = FALSE)
