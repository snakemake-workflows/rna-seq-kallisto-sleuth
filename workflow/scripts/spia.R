suppressPackageStartupMessages({
  library("httr") # only needed for the ssl_verifypeer = FALSE hack below
  library("SPIA")
  library("graphite")
})

# provides library("tidyverse") and function get_prefix_col()
# the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source( file.path(snakemake@scriptdir, 'common.R') )

options(Ncpus = snakemake@threads)

# TODO: This creates a security vulnerability, but it seems like
# something is currently wrong with the certificate chain for the
# downloads that bioconductor-graphite performs. To reproduce:
# `curl -vs "https://graphiteweb.bio.unipd.it/pathways/12/mmusculus.rda"`
httr::set_config(config(ssl_verifypeer = FALSE))

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
