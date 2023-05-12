log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")
library("tidyverse")
library("fs")
library("grid")
library("gridExtra")
library("dplyr")
model <- snakemake@params[["model"]]

sleuth_object <- sleuth_load(snakemake@input[[1]])

sleuth_object <- sleuth_lrt(sleuth_object, "reduced", "full")

# plot mean-variance
mean_var_plot <- plot_mean_var(sleuth_object, 
                               which_model = "full", 
                               point_alpha = 0.4,
                               point_size = 2, 
                               point_colors = c("black", "dodgerblue"),
                               smooth_alpha = 1, 
                               smooth_size = 0.75, 
                               smooth_color = "red")
ggsave(snakemake@output[["mean_var_plot"]], mean_var_plot)

write_results <- function(so, mode, output, output_all) {
    so$pval_aggregate <- FALSE
    if(mode == "aggregate") {
      # workaround the following bug-request:
      # https://github.com/pachterlab/sleuth/pull/240
      # TODO renaming can be removed when fixed
      g_col <- so$gene_column
      so$gene_column <- NULL
      so$pval_aggregate <- TRUE
      so$gene_column <- g_col
    }

    plot_model <- snakemake@wildcards[["model"]]

    # list for qq-plots to make a multipage pdf-file as output
    qq_list <- list()

    print("Performing likelihood ratio test")
    all <- sleuth_results(so, "reduced:full", "lrt", show_all = TRUE, pval_aggregate = so$pval_aggregate) %>%
            arrange(pval)
    all<- all %>% mutate_if(is.numeric, signif, 3)
    covariates <- colnames(design_matrix(so, "full"))
    covariates <- covariates[covariates != "(Intercept)"]

    # iterate over all covariates and perform wald test in order to obtain beta estimates
    if(!so$pval_aggregate) {

      # lists for volcano and ma-plots to make a multipage pdf-file as output
      volcano_list <- list()
      ma_list <- list()

      for(covariate in covariates) {
            print(str_c("Performing wald test for ", covariate))
            so <- sleuth_wt(so, covariate, "full")

            volc_plot_title <- str_c(plot_model, ": volcano plot for ", covariate)
            ma_plot_title <- str_c(plot_model, ": ma-plot for ", covariate)
            qq_plot_title <- str_c(plot_model, ": qq-plot from wald test for ", covariate)

            # volcano plot
            print(str_c("Performing volcano plot for ", covariate))
            volcano <- plot_volcano(so, covariate, "wt", "full",
                                    sig_level = snakemake@params[["sig_level_volcano"]], 
                                    point_alpha = 0.2, 
                                    sig_color = "red",
                                    highlight = NULL) +
              ggtitle(volc_plot_title) +
              theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
            volcano_list[[volc_plot_title]] <- volcano

            # ma-plot
            print(str_c("Performing ma-plot for ", covariate))
            ma_plot <- plot_ma(so, covariate, "wt", "full",
                               sig_level = snakemake@params[["sig_level_ma"]], 
                               point_alpha = 0.2, 
                               sig_color = "red",
                               highlight = NULL, 
                               highlight_color = "green") +
              ggtitle(ma_plot_title) +
              theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
            ma_list[[ma_plot_title]] <- ma_plot

            # qq-plots from wald test
            print(str_c("Performing qq-plot from wald test for ", covariate))
            qq_plot <- plot_qq(so, covariate, "wt", "full", 
                               sig_level = snakemake@params[["sig_level_qq"]],
                               point_alpha = 0.2, 
                               sig_color = "red", 
                               highlight = NULL, 
                               highlight_color = "green",
                               line_color = "blue") +
              ggtitle(qq_plot_title) +
              theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
            qq_list[[qq_plot_title]] <- qq_plot


            beta_col_name <- str_c("b", covariate, sep = "_")
            beta_se_col_name <- str_c(beta_col_name, "se", sep = "_")
            qval_name <- str_c("qval", covariate, sep = "_")
            all_wald <- sleuth_results(so, covariate, "wt", show_all = TRUE, pval_aggregate = FALSE) %>%
                        dplyr::select( target_id = target_id,
                                !!qval_name := qval,
                                !!beta_col_name := b,
                                !!beta_se_col_name := se_b)
            signed_pi_col_name <- str_c("signed_pi_value", covariate, sep = "_")
            all <- inner_join(all, all_wald, by = "target_id") %>%
	              # calculate a signed version of the pi-value from:
                    # https://dx.doi.org/10.1093/bioinformatics/btr671
	              # e.g. useful for GSEA ranking
	              mutate( !!signed_pi_col_name := -log10(pval) * !!sym(beta_col_name) )
      }
      # saving volcano plots
      marrange_volcano <- marrangeGrob(grobs=volcano_list, nrow=1, ncol=1, top = NULL)
      ggsave(snakemake@output[["volcano_plots"]], plot = marrange_volcano, width = 14)

      # saving ma-plots
      marrange_ma <- marrangeGrob(grobs=ma_list, nrow=1, ncol=1, top = NULL)
      ggsave(snakemake@output[["ma_plots"]], plot = marrange_ma, width = 14)
    }

    if(mode == "mostsignificant") {
        # for each gene, select the most significant transcript (this is equivalent to sleuth_gene_table, but with bug fixes)
      all <- all %>%
                drop_na(ens_gene) %>%
                group_by(ens_gene) %>%
                filter( qval == min(qval, na.rm = TRUE) ) %>%
                # ties in qval (e.g. min(qval) == 1) can lead to multiple entries per gene
                filter( pval == min(pval, na.rm = TRUE) ) %>%
                # for min(qval) == 1, then usually also min(pval) == 1, and
                # we need something else to select a unique entry per gene
                filter( mean_obs == max(mean_obs, na.rm = TRUE) ) %>%
                # sometimes we still get two transcript variants with exactly
                # the same values, i.e. they have exactly the same sequence
                # but (slightly) different annotations -- then we retain a string
                # with a comma-separated list of them
                mutate(target_id = str_c(target_id, collapse=",")) %>%
                distinct() %>%
                # useful sort for scrolling through output by increasing q-values
                arrange(qval)
      all <- all %>% mutate_if(is.numeric, signif, 3)
      # qq-plot from likelihood ratio test
      print(str_c("Performing qq-plot from likelihood ratio test"))
      qq_plot_title_trans <- str_c(plot_model, ": qq-plot from likelihood ratio test")
      qq_plot_trans <- plot_qq(so, 
                                test = 'reduced:full', 
                                test_type = 'lrt', 
                                sig_level = snakemake@params[["sig_level_qq"]],
                                point_alpha = 0.2, 
                                sig_color = "red", 
                                highlight = NULL, 
                                highlight_color = "green",
                                line_color = "blue") +
        ggtitle(qq_plot_title_trans) +
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
        qq_list[[qq_plot_title_trans]] <- qq_plot_trans
    } else if (mode == "canonical") {
      all <- all %>%
                drop_na(canonical) %>%
                filter(canonical)
      if (nrow(all) == 0) {
        stop("No canonical transcripts found (does ensembl support canonical transcript annotation for your species?")
      }
      # Control FDR again, because we have less tests now.
      all$qval <- p.adjust(all$pval, method = "BH")
      all <- all %>% mutate_if(is.numeric, signif, 3)
    } else if (mode == "custom") {
      # load custom ID list
      id_version_pattern <- "\\.\\d+$"
      ids <- read_tsv(snakemake@params[["representative_transcripts"]], col_names = "ID")$ID
      ids <- str_replace(ids, id_version_pattern, "")
      all <- all %>% 
        mutate(target_id_stem = str_replace(target_id, id_version_pattern, "")) %>% 
        filter(target_id_stem %in% ids) %>% 
        mutate(target_id_stem = NULL)
      all <- all %>% mutate_if(is.numeric, signif, 3)
      if (nrow(all) == 0) {
        stop("The given list of representative transcript ids does not match any of the transcript ids of the chosen species.")
      }
    }

    # saving qq-plots
    marrange_qq <- marrangeGrob(grobs=qq_list, nrow=1, ncol=1, top = NULL)
    ggsave(snakemake@output[["qq_plots"]], plot = marrange_qq, width = 14)

    write_rds(all, path = output_all, compress = "none")
    # add sample expressions
    all <- all %>% left_join(as_tibble(sleuth_to_matrix(so, "obs_norm", "est_counts"), rownames="target_id"))
    all <- all %>% mutate_if(is.numeric, signif, 3)
    write_tsv(all, path = output, escape = "none")
}

write_results(sleuth_object, "transcripts", snakemake@output[["transcripts"]], snakemake@output[["transcripts_rds"]])
write_results(sleuth_object, "aggregate", snakemake@output[["genes_aggregated"]], snakemake@output[["genes_aggregated_rds"]])

repr_trans <- snakemake@params[["representative_transcripts"]]
if (repr_trans == "canonical") {
  write_results(sleuth_object, "canonical", snakemake@output[["genes_representative"]], snakemake@output[["genes_representative_rds"]])
} else if (repr_trans == "mostsignificant") {
  write_results(sleuth_object, "mostsignificant", snakemake@output[["genes_representative"]], snakemake@output[["genes_representative_rds"]])
} else {
  write_results(sleuth_object, "custom", snakemake@output[["genes_representative"]], snakemake@output[["genes_representative_rds"]])
}

