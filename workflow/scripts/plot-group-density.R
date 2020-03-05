log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")


so <- sleuth_load(snakemake@input[[1]])

plot_group_density(so,
                   use_filtered = TRUE,
                   units = "est_counts",
                   trans = "log",
                   grouping = setdiff(colnames(so$sample_to_covariates), "sample"),
                   offset = 1)
ggsave(snakemake@output[[1]])
