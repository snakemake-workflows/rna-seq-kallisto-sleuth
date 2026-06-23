log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(sleuth)
library(tidyverse)

  # - bioconductor-xcell2  zwei schritte 1 sleuth dann decon um parallel 
  # - bioconductor-complexheatmap
  # - r-seriation

so <- sleuth_load(snakemake@input[[1]])

# transcript-level, sleuth-normalized, pre-log counts
# rows = target_id (transcript), columns = sample
norm_counts <- sleuth_to_matrix(so, "obs_norm", "est_counts")

# map transcripts -> genes using the sleuth object's own target_mapping,
target_mapping <- so$target_mapping %>%
  dplyr::select(target_id, ens_gene) %>%
  drop_na(ens_gene)

gene_counts <- as.data.frame(norm_counts) %>%
  rownames_to_column(var = "target_id") %>%
  inner_join(target_mapping, by = "target_id") %>%
  dplyr::select(-target_id) %>%
  group_by(ens_gene) %>%
  summarise(across(everything(), sum), .groups = "drop")

# drop genes that are all-zero across every sample
sample_cols <- setdiff(colnames(gene_counts), "ens_gene")
nonzero_mask <- rowSums(gene_counts[, sample_cols, drop = FALSE]) > 0
gene_counts <- gene_counts[nonzero_mask, ]

# log2(count + 3) transform, as requested for downstream analysis 
# article verweis 
gene_counts[, sample_cols] <- log2(gene_counts[, sample_cols] + 3)

write_tsv(gene_counts, snakemake@output[["gene_counts"]])

### now I have to get the gene names used in the ref sets 
### here the ref set can be down loaded 
### https://github.com/AlmogAngel/xCell2/tree/master/data