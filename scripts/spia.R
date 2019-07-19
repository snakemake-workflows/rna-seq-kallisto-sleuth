library(sleuth)
library(SPIA)
library(graphite)
library(AnnotationDbi)
library(tidyverse)

options(Ncpus = snakemake@threads)

covariate <- snakemake@params[["covariate"]]

prefix_ens_ids <- function(ids) {
    paste("ENSEMBL", ids, sep = ":")
}

db <- pathways(snakemake@params[["species"]], "reactome")
db <- convertIdentifiers(db, "ENSEMBL")

prepareSPIA(db, "reactome")


diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>% drop_na(ens_gene) %>% mutate(ens_gene = str_c("ENSEMBL", ens_gene))
universe <- diffexp %>% pull(ens_gene)
diffexp <- diffexp %>% filter(qval <= 0.05)

# get logFC equivalent (the sum of beta scores of covariates of interest)
beta_col <- paste("b", covariate, sep = "_")

cols <- colnames(diffexp)
suffixes <- c("", "1", ".0")
found <- FALSE
for(suffix in suffixes) {
    beta_col <- paste0(beta_col, suffix)
    if(beta_col %in% cols) {
        found <- TRUE
        break
    }
}
if(!found) {
    stop(paste0("Invalid covariate ", covariate, ", not found in diffexp table."))
}

beta <- diffexp %>% select(ens_gene, beta_col) %>% deframe()

res <- runSPIA(de = beta, all = universe, "reactome")

write_tsv(res, snakemake@output[[1]])
