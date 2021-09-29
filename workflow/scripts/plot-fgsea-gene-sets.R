log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("fgsea")

# provides library("tidyverse") and function get_prefix_col()
# the latter requires snakemake@output[["samples"]] and
# snakemake@params[["covariate"]]
source(snakemake@params[["common_src"]])


gene_sets <- gmtPathways(snakemake@input[["gene_sets"]])
sig_gene_sets <- read_tsv(snakemake@input[["sig_gene_sets"]])
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

dir.create( snakemake@output[[1]] )

for ( set in (sig_gene_sets %>% pull(pathway)) ) {
  # plot gene set enrichment
  if (length(gene_sets[[set]]) == 0 || is.na(ranked_genes[as.vector(gene_sets[[set]])])) {
    next
  }
  p <- plotEnrichment(gene_sets[[set]], ranked_genes) +
         ggtitle(str_c("gene set: ", set),
           subtitle = str_c(
             sig_gene_sets %>% filter(pathway == set) %>% pull(size)," genes;  ",
             "p-value (BH-adjusted): ", sig_gene_sets %>% filter(pathway ==  set) %>% pull(padj), "\n",
             "normalized enrichment score (NES):", sig_gene_sets %>% filter(pathway == set) %>% pull(NES)
             )
         ) +
         xlab("gene rank") +
         theme_bw( base_size = 16 )
  setname <- gsub("[-%/:,'\\. ]", "_", set)
  fname <- str_c(snakemake@wildcards[["model"]], ".", setname, ".gene-set-plot.pdf")
  ggsave(
    file = file.path( snakemake@output[[1]], fname ),
    width = 10,
    height = 7
  )
}
