log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("fgsea")
library("AnnotationDbi")

# provides library("tidyverse") and function get_prefix_col()
# the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source( file.path(snakemake@scriptdir, 'common.R') )

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


gene_sets <- gmtPathways(snakemake@input[["gene_sets"]])
diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
                  drop_na(ext_gene) %>%
                  mutate(ext_gene = str_to_upper(ext_gene)) %>%
                  group_by(ext_gene) %>%
                	filter( qval == min(qval, na.rm = TRUE) ) %>%
                	mutate(target_id = str_c(target_id, collapse=",")) %>%
                	mutate(ens_gene = str_c(ens_gene, collapse=",")) %>%
                  distinct()

signed_pi <- get_prefix_col("signed_pi_value", colnames(diffexp))

ranked_genes <- diffexp %>%
                  dplyr::select(ext_gene, !!signed_pi) %>%
                  deframe()

# get and write out rank values that are tied -- a way to check up on respecitve warnings
rank_ties <- enframe(ranked_genes) %>%
               mutate(dup = duplicated(value) | duplicated(value, fromLast = TRUE) ) %>%
               filter(dup == TRUE) %>%
               dplyr::select(ext_gene = name, !!signed_pi := value)
write_tsv(rank_ties, snakemake@output[["rank_ties"]])

fgsea_res <- fgsea(pathways = gene_sets, 
                    stats = ranked_genes,
                    minSize=10,
                    maxSize=700,
                    nproc=snakemake@threads,
                    nperm=snakemake@params[["nperm"]]
                    ) %>%
                as_tibble() 

# Annotation is impossible without any entries, so then just write out empty files
if ( (fgsea_res %>% count() %>% pull(n)) == 0 ) {

    write_tsv(fgsea_res, path = snakemake@output[["enrichment"]])
    write_tsv(fgsea_res, path = snakemake@output[["significant"]])

} else {

    # add further annotation
    annotated <- fgsea_res %>%
                    unnest(leadingEdge) %>%
                    mutate(
                        leading_edge_symbol = str_to_title(leadingEdge),
                        leading_edge_entrez_id = mapIds(leading_edge_symbol, x=get(pkg), keytype="SYMBOL", column="ENTREZID"),
                        leading_edge_ens_gene = mapIds(leading_edge_symbol, x=get(pkg), keytype="SYMBOL", column="ENSEMBL")
                        ) %>%
                    group_by(pathway) %>%
                    summarise(
                        leadingEdge = str_c(leadingEdge, collapse = ','),
                        leading_edge_symbol = str_c(leading_edge_symbol, collapse = ','),
                        leading_edge_entrez_id = str_c(leading_edge_entrez_id, collapse = ','),
                        leading_edge_ens_gene = str_c(leading_edge_ens_gene, collapse = ',')
                    ) %>%
                    inner_join(fgsea_res %>% select(-leadingEdge), by = "pathway") %>%
                    select(-leadingEdge, -leading_edge_symbol,
                           -leading_edge_entrez_id, -leading_edge_ens_gene,
                            leading_edge_symbol, leading_edge_ens_gene,
                            leading_edge_entrez_id, leadingEdge)
    
    # write out fgsea results for all gene sets
    write_tsv(annotated, path = snakemake@output[["enrichment"]])
    
    # select significant pathways
    sig_gene_sets <- annotated %>%
                       filter( padj < snakemake@params[["gene_set_fdr"]] )
    
    # write out fgsea results for gene sets found to be significant
    write_tsv(sig_gene_sets, path = snakemake@output[["significant"]])
}

height = .7 * (length(gene_sets) + 2)

# table plot of all gene sets
tg <- plotGseaTable(
            pathway = gene_sets,
            stats = ranked_genes,
            fgseaRes = fgsea_res,
            gseaParam = 1,
            render = FALSE
        )
ggsave(filename = snakemake@output[["plot"]], plot = tg, width = 12, height = height, limitsize=FALSE)
