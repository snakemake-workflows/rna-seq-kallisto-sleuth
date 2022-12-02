log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("fgsea")
library(snakemake@params[["bioc_species_pkg"]], character.only = TRUE)

# provides library("tidyverse") and functions load_bioconductor_package() and
# get_prefix_col(), the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source(snakemake@params[["common_src"]])

gene_sets <- gmtPathways(snakemake@input[["gene_sets"]])
diffexp <- read_tsv(snakemake@input[["diffexp"]]) %>%
                  drop_na(ext_gene) %>%
                  mutate(ext_gene = str_to_upper(ext_gene)) %>%
                  # resolve multiple occurences of the same ext_gene, usually
                  # meaning that several ENSEMBLE genes map to the same gene
                  # symbol -- so we may loose resolution here, as long as gene
                  # symbols are used
                  group_by(ext_gene) %>%
                    filter( qval == min(qval, na.rm = TRUE) ) %>%
                    filter( pval == min(pval, na.rm = TRUE) ) %>%
                    # for the case of min(qval) == 1 and min(pval) == 1, we
                    # need something else to select a unique entry per gene
                    filter( mean_obs == max(mean_obs, na.rm = TRUE) ) %>%
                    mutate(target_id = str_c(target_id, collapse=",")) %>%
                    mutate(ens_gene = str_c(ens_gene, collapse=",")) %>%
                  distinct()

signed_pi <- get_prefix_col("signed_pi_value", colnames(diffexp))

ranked_genes <- diffexp %>%
                  dplyr::select(ext_gene, !!signed_pi) %>%
                  deframe()

# get and write out rank values that are tied -- a way to check up on respective warnings
rank_ties <- enframe(ranked_genes) %>%
               mutate(dup = duplicated(value) | duplicated(value, fromLast = TRUE) ) %>%
               filter(dup == TRUE) %>%
               dplyr::select(ext_gene = name, !!signed_pi := value)
write_tsv(rank_ties, snakemake@output[["rank_ties"]])

print(gene_sets)
print(ranked_genes)

fgsea_res <- fgsea(pathways = gene_sets,
                    stats = ranked_genes,
                    minSize=10,
                    maxSize=700,
                    nproc=snakemake@threads,
                    eps=snakemake@params[["eps"]]
                    ) %>%
                as_tibble()

# Annotation is impossible without any entries, so then just write out empty files
if ( (fgsea_res %>% count() %>% pull(n)) == 0 ) {

    # leadingEdge cannot be a list, and the empty output should at least
    # have the correct columns in the correct order
    unnested <- fgsea_res %>%
                    unnest(leadingEdge) %>%
                    add_column(
                        leading_edge_symbol = NA,
                        leading_edge_ens_gene = NA,
                        leading_edge_entrez_id = NA
                    ) %>%
                    dplyr::relocate(
                        c(
                          leading_edge_symbol,
                          leading_edge_ens_gene,
                          leading_edge_entrez_id,
                          leadingEdge
                        ),
                        .after = last_col()
                    )

    write_tsv(unnested, file = snakemake@output[["enrichment"]])
    write_tsv(unnested, file = snakemake@output[["significant"]])

} else {

    # add further annotation
    annotated <- fgsea_res %>%
                    unnest(leadingEdge) %>%
                    mutate(
                        leading_edge_symbol = str_to_title(leadingEdge),
                        leading_edge_entrez_id = mapIds(leading_edge_symbol, x=get(snakemake@params[["bioc_species_pkg"]]), keytype="SYMBOL", column="ENTREZID"),
                        leading_edge_ens_gene = mapIds(leading_edge_symbol, x=get(snakemake@params[["bioc_species_pkg"]]), keytype="SYMBOL", column="ENSEMBL")
                        ) %>%
                    group_by(pathway) %>%
                    summarise(
                        leadingEdge = str_c(leadingEdge, collapse = ','),
                        leading_edge_symbol = str_c(leading_edge_symbol, collapse = ','),
                        leading_edge_entrez_id = str_c(leading_edge_entrez_id, collapse = ','),
                        leading_edge_ens_gene = str_c(leading_edge_ens_gene, collapse = ',')
                    ) %>%
                    inner_join(fgsea_res %>% dplyr::select(-leadingEdge), by = "pathway") %>%
                    dplyr::relocate(
                        c(
                          leading_edge_symbol,
                          leading_edge_ens_gene,
                          leading_edge_entrez_id,
                          leadingEdge
                        ),
                        .after = last_col()
                    )

    # write out fgsea results for all gene sets
    write_tsv(annotated, file = snakemake@output[["enrichment"]])

    # select significant pathways
    sig_gene_sets <- annotated %>%
                       filter( padj < snakemake@params[["gene_set_fdr"]] )

    # write out fgsea results for gene sets found to be significant
    write_tsv(sig_gene_sets, file = snakemake@output[["significant"]])
}

# select significant pathways
top_pathways <- fgsea_res %>% arrange(padj) %>% head(n=1000) %>% filter(padj <= snakemake@params[["gene_set_fdr"]]) %>% arrange(-NES) %>% pull(pathway)
selected_gene_sets <- gene_sets[top_pathways]

height = .7 * (length(selected_gene_sets) + 2)

# table plot of all gene sets
tg <- plotGseaTable(
            pathways = selected_gene_sets,
            stats = ranked_genes,
            fgseaRes = fgsea_res,
            gseaParam = 1,
            render = FALSE
        )
ggsave(filename = snakemake@output[["plot"]], plot = tg, width = 15, height = height, limitsize=FALSE)

collapsed_pathways <- collapsePathways(fgsea_res %>% arrange(pval) %>% filter(padj <= snakemake@params[["gene_set_fdr"]]), gene_sets, ranked_genes)
main_pathways <- fgsea_res %>% filter(pathway %in% collapsed_pathways$mainPathways) %>% arrange(-NES) %>% pull(pathway)
selected_gene_sets <- gene_sets[main_pathways]
height = .7 * (length(selected_gene_sets) + 2)

# table plot of all gene sets
tg <- plotGseaTable(
            pathways = selected_gene_sets,
            stats = ranked_genes,
            fgseaRes = fgsea_res,
            gseaParam = 1,
            render = FALSE
        )
ggsave(filename = snakemake@output[["plot_collapsed"]], plot = tg, width = 15, height = height, limitsize=FALSE)
