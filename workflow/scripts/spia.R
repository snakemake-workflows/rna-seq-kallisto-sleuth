log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("SPIA")
library("graphite")
library(snakemake@params[["bioc_species_pkg"]], character.only = TRUE)

# provides library("tidyverse") and get_prefix_col(), where the latter requires
# snakemake@input[["samples"]] and snakemake@params[["covariate"]]
source(snakemake@input[["common_src"]])

pw_db <- snakemake@params[["pathway_db"]]
db <- readRDS(snakemake@input[["spia_db"]])

orgDb_name <- snakemake@params[["orgDb"]]
orgDb <- if (orgDb_name == "NA") NA else get(orgDb_name)
print(typeof(orgDb))

options(Ncpus = snakemake@threads)

diffexp <- read_tsv(snakemake@input[["diffexp"]]) |>
  drop_na(ens_gene) |>
  mutate(ens_gene = str_c("ENSEMBL:", ens_gene))

universe <- diffexp |>
  dplyr::select(ens_gene) |>
  distinct() |>
  pull(ens_gene)

sig_genes <- diffexp |>
  filter(qval <= 0.05)

columns <- c(
  "Name",
  "number of genes on the pathway",
  "number of DE genes per pathway",
  "p-value for at least NDE genes",
  "total perturbation accumulation",
  "p-value to observe a total accumulation",
  "Combined p-value",
  "Combined FDR",
  "Combined Bonferroni p-values",
  "Status",
  "pathway id",
  "gene_ratio",
  "study_items_sig_terms"
)

if (nrow(sig_genes) == 0) {
  # the best hack for an empty tibble from a column specification I could find
  res <- read_csv("\n", col_names = columns)
  write_tsv(res, snakemake@output[["table"]])
  # Create empty plot
  pdf(snakemake@output[["plots"]])
  plot.new()
  text(0.5, 0.5, "No significant genes found", cex = 1.5)
  dev.off()
} else {
  # get logFC equivalent (the sum of beta scores of covariates of interest)

  beta_col <- get_prefix_col("b", colnames(sig_genes))

  beta <- sig_genes |>
    dplyr::select(ens_gene, !!beta_col) |>
    deframe()

  t <- tempdir(check = TRUE)
  olddir <- getwd()
  setwd(t)
  prepareSPIA(db, pw_db)
  res <- runSPIA(
    de = beta,
    all = universe,
    pw_db,
    plots = TRUE,
    verbose = TRUE
  )
  setwd(olddir)

  file.copy(
    file.path(t, "SPIAPerturbationPlots.pdf"),
    snakemake@output[["plots"]]
  )
  pathway_names <- db[res$Name]
  if (length(pathway_names) > 0) {
    pathway_ids_tibble <- pathway_names@entries |>
      map(slot, "id") |>
      unlist() |>
      as_tibble(
        rownames="pathway_name"
      ) |>
      rename(
        `pathway id` = value
      )

    extract_genes_from_pathway <- function(pathway) {
      # Extract genes from protEdges and mixedEdges
      genes_protEdges <- unique(c(pathway@protEdges$src, pathway@protEdges$dest))
      genes_mixedEdges <- unique(c(pathway@mixedEdges$src[pathway@mixedEdges$src_type == "ENSEMBL"], 
                                  pathway@mixedEdges$dest[pathway@mixedEdges$dest_type == "ENSEMBL"]))
      
      # Combine all gene IDs
      all_genes <- unique(c(genes_protEdges, genes_mixedEdges))
      
      # Map ENSMUSB values to gene names
      external_gene_names <- mapIds(orgDb, 
                                    keys = all_genes, 
                                    column = "SYMBOL", 
                                    keytype = "ENSEMBL",
                                    multiVals = "first")
      
      # Check if external_gene_names is empty or NULL
      if (is.null(external_gene_names) || length(external_gene_names) == 0) {
          return("")
      }

      # Filter diffexp for relevant gene names
      filtered_diffexp <- diffexp %>%
          filter(ext_gene %in% unname(external_gene_names)) %>%
          select(ext_gene, beta_val = !!beta_col) %>%
          drop_na()
      
      # Check if filtered_diffexp is empty
      if (nrow(filtered_diffexp) == 0) {
          return("")
      }

      # Create result string
      result_string <- filtered_diffexp %>%
          mutate(pair = paste(ext_gene, beta_val, sep = ":")) %>%
          pull(pair) %>%
          paste(collapse = ", ")
      
      return(result_string)
    }

    extract_genes_from_pathway_v <- Vectorize(extract_genes_from_pathway)

    # Create new column with genes in pathway together with their beta_vals
    if (!is.na(orgDb)) {
      res <- res %>%
          mutate(study_items_sig_terms = extract_genes_from_pathway_v(lapply(res$Name, function(name) db[[name]])))
    }
    else {
      res[, 'study_items_sig_terms'] = NA

    }

    final_res <- as_tibble(res) |>
      left_join(
        pathway_ids_tibble,
        join_by(Name == pathway_name)
      ) |>
      rename(
        "number of genes on the pathway" = "pSize",
        "number of DE genes per pathway" = "NDE",
        "p-value for at least NDE genes" = "pNDE",
        "total perturbation accumulation" = "tA",
        "p-value to observe a total accumulation" = "pPERT",
        "Combined p-value" = "pG",
        "Combined FDR" = "pGFdr",
        "Combined Bonferroni p-values" = "pGFWER"
      ) |>
      mutate(
        gene_ratio = str_c("(", `number of DE genes per pathway`, ", ", `number of genes on the pathway`, ")")
      ) |>
      dplyr::select(all_of(columns)
      ) |>
      arrange(`Combined FDR`)
     
    write_tsv(final_res, snakemake@output[["table"]])
  } else {
    # the best hack for an empty tibble from a column specification I could find
    emtpy_data_frame <- read_csv("\n", col_names = columns)
    write_tsv(emtpy_data_frame, snakemake@output[["table"]])
  }
}
