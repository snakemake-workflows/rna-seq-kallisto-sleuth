log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("limma")
library("tidyverse")

so <- sleuth_load(snakemake@input[[1]])

# get normed counts
norm_counts <- sleuth_to_matrix(so, "obs_norm", "est_counts")

# log transform
log_norm_counts <- log2(norm_counts + 1)

# reorder columns to match covariates
log_norm_counts <- log_norm_counts[, so$sample_to_covariates$sample, drop=FALSE]

# obtain covariates
model <- snakemake@params[["model"]]
model_nobatch <- paste("~", model[["primary_variable"]], sep="")

covariates <- model.matrix(as.formula(model[["reduced"]]), data = so$sample_to_covariates)
design <- model.matrix(as.formula(model_nobatch), data = so$sample_to_covariates)

stopifnot(so$sample_to_covariates$sample == colnames(log_norm_counts))

# remove batch effects
final_counts <- removeBatchEffect(log_norm_counts, covariates = covariates, design = design)

target_mapping <- so$target_mapping
rownames(target_mapping) <- target_mapping$target_id

final_counts <- rownames_to_column(as.data.frame(final_counts), var="transcript") %>% 
                add_column(
                    gene = target_mapping[rownames(final_counts), "ext_gene"],
                    .after = "transcript"
                )


write_tsv(final_counts, snakemake@output[[1]])

