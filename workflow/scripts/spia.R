suppressPackageStartupMessages({
  library("SPIA")
  library("graphite")
  library("AnnotationDbi")
})

# provides library("tidyverse") and function get_prefix_col()
# the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source( file.path(snakemake@scriptdir, 'common.R') )

options(Ncpus = snakemake@threads)

db <- pathways(snakemake@params[["species"]], "reactome")
db <- convertIdentifiers(db, "ENSEMBL")

prepareSPIA(db, "reactome")


diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
            drop_na(ens_gene) %>%
            mutate(ens_gene = str_c("ENSEMBL", ens_gene))
universe <- diffexp %>% pull(var = ens_gene)
sig_genes <- diffexp %>% filter(qval <= 0.05)

# get logFC equivalent (the sum of beta scores of covariates of interest)

beta_col <- get_prefix_col("b", colnames(sig_genes))

beta <- sig_genes %>%
            select(ens_gene, !!beta_col) %>%
            deframe()

res <- runSPIA(de = beta, all = universe, "reactome", plots = TRUE)

write_tsv(res, snakemake@output[[1]])
