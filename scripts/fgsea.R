suppressPackageStartupMessages({
  library("fgsea")
  library("AnnotationDbi")
})

# provides library("tidyverse") and function get_beta_col()
source('scripts/common.R')

covariate <- snakemake@params[["covariate"]]

# get the species for adding annotations below
species <- str_to_title( str_sub(snakemake@params[["species"]], 1, 2) )

# construct the respective package name
pkg <- str_c("org.", species, ".eg.db")

# check if the package is installed in the fgsea environment and load,
# give a useful error message otherwise; installed per default are:
# * org.Mm.eg.db
# * org.Hs.eg.db
if ( pkg %in% rownames( installed.packages() ) ) {
    library(pkg, character.only = TRUE)
} else {
    stop(
        str_c(
            "\n",
            "Package not installed: ", pkg, "\n",
            "Package name inferred from species name: ", species, "\n",
            "Check if package 'bioconductor-", pkg, "' exists and add to envs/fgsea.yaml\n",
            "\n"
        )
    )
}


# load gene set specified via wildcard
gene_sets <- gmtPathways(snakemake@input[["gene_sets"]])
diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
                  drop_na(ext_gene) %>%
                  mutate(ext_gene = str_to_upper(ext_gene)) %>%
			group_by(ext_gene) %>%
                	filter( qval == min(qval, na.rm = TRUE) ) %>%
                	mutate(target_id = str_c(target_id, collapse=",")) %>%
                	mutate(ens_gene = str_c(ens_gene, collapse=",")) %>%
			distinct()

beta_col <- get_beta_col(covariate, colnames(diffexp))

ranked_genes <- diffexp %>%
                  dplyr::select(ext_gene, !!beta_col) %>%
                  deframe()

# get and write out rank values that are tied -- a way to check up on respecitve warnings
rank_ties <- enframe(ranked_genes) %>%
               mutate(dup = duplicated(value) | duplicated(value, fromLast = TRUE) ) %>%
               filter(dup == TRUE) %>%
               dplyr::select(ext_gene = name, !!beta_col := value)
write_tsv(rank_ties, snakemake@output[["rank_ties"]])

fgsea_res <- fgsea(pathways = gene_sets, 
                    stats = ranked_genes,
                    minSize=10,
                    maxSize=700,
                    nproc=snakemake@threads,
                    nperm=strtoi(snakemake@wildcards[["nperm"]])
                ) %>%
                as_tibble() %>%
                # add further annotation
                mutate(
                    symbols = sapply(leadingEdge, str_to_title),
                    entrez_ids = str_c(
                        lapply(symbols, mapIds, x=get(pkg), keytype="SYMBOL", column="ENTREZID"),
                        sep = ','),
                    ens_gene = str_c(
                        lapply(symbols, mapIds, x=get(pkg), keytype="SYMBOL", column="ENSEMBL"),
                        sep = ','),
                    symbols = str_c(symbols, sep = ','),
                    leadingEdge = str_c(leadingEdge, sep = ',')
                )

# plot results for all pathways
write_tsv(fgsea_res, path = snakemake@output[["enrichment"]])

# select significant pathways
sig_gene_sets <- fgsea_res %>%
                   filter( padj < snakemake@params[["gene_set_fdr"]] )

# plot significant pathways
dir.create(snakemake@output[["plot"]])
for (set in sig_gene_sets %>% pull(pathway)) {
    p <- plotEnrichment(gene_sets[[set]], ranked_genes) +
            ggtitle(str_c("gene set: ", set),
                subtitle = str_c(
                    sig_gene_sets %>% filter(pathway ==  set) %>% pull(size), "genes;  ",
                    "p-value (BH-adjusted): ", sig_gene_sets %>% filter(pathway ==  set) %>% pull(padj), "\n",
                    "normalized enrichment score (NES):", sig_gene_sets %>% filter(pathway == set) %>% pull(NES)
                    )
            ) +
            xlab("gene rank") + 
            theme_bw(base_size=16)
    ggsave(file.path(snakemake@output[["plot"]], str_c("fgsea_enrichment-plot_significant-gene-set_", set, ".pdf")), width=10, height=7)
}

if (length(sig_gene_sets) == 0) {
    p <- ggplot() +
            ggtitle("No gene sets found significant after\nBenjamini-Hochberg false discovery rate control.")
    ggsave(file.path(snakemake@output[["plot"]], "dummy.pdf"), width=10, height=7)
}

