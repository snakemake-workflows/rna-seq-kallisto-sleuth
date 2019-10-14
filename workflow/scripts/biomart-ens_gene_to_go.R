log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("biomaRt")
library("tidyverse")

# create an ensembl biomart for the species specified in the params field
mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
  host = 'ensembl.org')

# get all ensembl_gene_id - go_id pairs and collapse into key->values mapping:
#     ensembl_gene_id -> go_id;go_id;go_id
ens_gene_to_go <- biomaRt::getBM(attributes = c("ensembl_gene_id", "go_id"), mart = mart) %>%
                    # rows with empty go_id are useless and will lead to extra semicolons below
                    filter(go_id != "") %>%
                    group_by(ensembl_gene_id) %>%
                    summarise(go_ids = str_c(go_id, collapse = ";"))

write_tsv(ens_gene_to_go, path = snakemake@output[[1]], col_names = FALSE)

