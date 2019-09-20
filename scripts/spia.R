library(sleuth)
library(SPIA)
library(graphite)
library(AnnotationDbi)
library(tidyverse)

options(Ncpus = snakemake@threads)

covariate <- snakemake@params[["covariate"]]

db <- pathways(snakemake@params[["species"]], "reactome")
db <- convertIdentifiers(db, "ENSEMBL")

prepareSPIA(db, "reactome")


diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
            drop_na(ens_gene) %>%
            mutate(ens_gene = str_c("ENSEMBL", ens_gene))
universe <- diffexp %>% pull(var = ens_gene)
sig_genes <- diffexp %>% filter(qval <= 0.05)

# get logFC equivalent (the sum of beta scores of covariates of interest)
beta_col <- str_c("b", covariate, sep = "_")

cols <- colnames(sig_genes)
suffixes <- c("", "1", ".0")
found <- FALSE
for(suffix in suffixes) {
    beta_col <- str_c(beta_col, suffix)
    if(beta_col %in% cols) {
        found <- TRUE
        break
    }
}
if(!found) {
    stop(str_c("Invalid covariate ", covariate, ", not found in diffexp table."))
}

beta <- sig_genes %>%
            select(ens_gene, !!beta_col) %>%
            deframe()

res <- runSPIA(de = beta, all = universe, "reactome", plots = TRUE)

write_tsv(res, snakemake@output[[1]])
