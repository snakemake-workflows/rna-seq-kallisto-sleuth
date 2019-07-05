library(SPIA)
library(graphite)

db <- pathways(snakemake@params[["species"]], "reactome")

prepareSPIA(db, "reactomeEx")

so <- sleuth_load(snakemake@input[["sleuth"]])

# get gene to transcript mapping
gene_table <- read.table(snakemake@input[["genes_to_transcripts"]], sep = "\t", header = TRUE)
all <- gene_table$ens_gene
significant <- gene_table[gene_table$qval <= 0.05, ]

# get differentially expressed genes
diffexp <- read.table(snakemake@input[["diffexp"]], sep = "\t", header = TRUE)
rownames(diffexp) <- diffexp$target_id
diffexp <- diffexp[significant$most_sig_transcript, ]

# obtain covariates of interest
covariates <- setdiff(colnames(design_matrix(so, "full")), colnames(design_matrix(so, "reduced")))

# get logFC equivalent (the sum of beta scores of covariates of interest)
beta <- rowSums(diffexp[, paste0("beta.", covariates)])
names(beta) <- diffexp$ens_gene

res <- runSPIA(de = beta, all = all, "reactomeEx")

write.table(res, file = snakemake@output[[1]], sep = "\t", col.names = NA, row.names = TRUE) 
