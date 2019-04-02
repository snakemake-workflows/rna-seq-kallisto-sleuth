suppressMessages({
  library("sleuth")
  library("biomaRt")
})

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = paste(snakemake@params[["species"]], "_gene_ensembl", sep = ""),
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

samples <- read.table(snakemake@input[["samples"]])
samples[, "path"] <- as.character(samples[, "path"])

so <- sleuth_prep(samples, extra_bootstrap_summary = TRUE, target_mapping = t2g, aggregation_column = "ens_gene")

sleuth_save(so, snakemake@output[[1]])
