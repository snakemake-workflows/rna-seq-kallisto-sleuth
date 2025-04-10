log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("limma")
library("tidyverse")

so <- sleuth_load(snakemake@input[[1]])

# get TPM values
tpm_matrix <- sleuth_to_matrix(so, "obs_norm", "tpm")

# reorder columns to match covariates
tpm_matrix <- tpm_matrix[, so$sample_to_covariates$sample, drop=FALSE]

target_mapping <- so$target_mapping
rownames(target_mapping) <- target_mapping$target_id

# add transcript and gene information
tpm_matrix <- rownames_to_column(as.data.frame(tpm_matrix), var="transcript") %>%
              add_column(
                gene = target_mapping[rownames(tpm_matrix), "ext_gene"],
                .after = "transcript"
              )

write_tsv(tpm_matrix, snakemake@output[[1]])
