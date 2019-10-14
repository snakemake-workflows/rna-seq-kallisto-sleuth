log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")

so <- sleuth_load(snakemake@input[[1]])

tpm <- sleuth_to_matrix(so, "obs_norm", "tpm")
target_mapping <- so$target_mapping
rownames(target_mapping) <- target_mapping$target_id
tpm <- cbind(ext_gene = target_mapping[rownames(tpm), "ext_gene"], tpm)

write.table(tpm, file = snakemake@output[[1]], col.names = NA, row.names = TRUE)

