library(sleuth)
library(SPIA)
library(graphite)
library(AnnotationDbi)

options(Ncpus = snakemake@threads)

covariate <- snakemake@params[["covariate"]]
#library(org.Hs.eg.db)
#print(columns(org.Hs.eg.db))
#annotation_db <- get(selectDb(snakemake@params[["species"]]))
#print(columns(annotation_db))

prefix_ens_ids <- function(ids) {
    paste("ENSEMBL", ids, sep = ":")
}

db <- pathways(snakemake@params[["species"]], "reactome")
db <- convertIdentifiers(db, "ENSEMBL")

prepareSPIA(db, "reactome")

so <- sleuth_load(snakemake@input[["sleuth"]])

diffexp <- read.table(snakemake@input[["diffexp"]], sep = "\t", header = TRUE)
diffexp <- diffexp[diffexp$qval <= 0.05, ]
diffexp$ens_gene <- prefix_ens_ids(diffexp$ens_gene)

#diffexp$entrez_gene <- mapIds(annotation_db, keys = diffexp$ens_gene, keytype = "ENSEMBLID", columns = "ENTREZID")

# get logFC equivalent (the sum of beta scores of covariates of interest)
beta_col <- paste("b", covariate, sep = "_")

cols <- colnames(diffexp)
suffixes <- c("", "1", ".0")
found <- FALSE
for(suffix in suffixes) {
    beta_col <- paste0(beta_col, suffix)
    if(beta_col %in% cols) {
        found <- TRUE
    }
}
if(!found) {
    stop(paste0("Invalid covariate ", covariate, ", not found in diffexp table."))
}

beta <- rowSums(diffexp[, beta_col, drop = FALSE])
names(beta) <- diffexp$ens_gene

res <- runSPIA(de = beta, all = prefix_ens_ids(unique(so$target_mapping$ens_gene)), "reactome")

write.table(res, file = snakemake@output[[1]], sep = "\t", col.names = NA, row.names = TRUE) 
