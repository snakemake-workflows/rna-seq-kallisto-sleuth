  library("tidyverse")
  library("ggpubr")
  library("IHW")
  
  # TODO: SET ABSOLUTE PATHS HERE
  ### this can be used as snakemake-input
  input_dir <- "..." # absolute path to directory of rna-seq-kallisto-sleuth
  
  ### this can be used as snakemake-output
  output_dir <- "..."
  
  output_dir_diagn_plots <- str_c(output_dir, "/diagnostic_plots") # directory for pdf's of diagnostic plots
  output_dir_result_plots <- str_c(output_dir, "/result_plots") # direcotry for pdf's of result plots
  
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } 
  if (!dir.exists(output_dir_diagn_plots)){
    dir.create(output_dir_diagn_plots)
  } 
  if (!dir.exists(output_dir_result_plots)){
    dir.create(output_dir_result_plots)
  } 
  
  setwd(input_dir)
  
  # this can be used as snakemake-param
  param_number.of.groups <- 7 
  
  # covariate: mean_obs in genes.aggregated
  genes.aggregated <- read_tsv(".test/results/tables/diffexp/model_X.genes-aggregated.diffexp.tsv")
  genes.aggregated <- na.omit(genes.aggregated) %>%
    mutate(adjusted_mean_obs = if_else(num_aggregated_transcripts > 0, (sum_mean_obs_counts / num_aggregated_transcripts), 0),
           grouping = groups_by_filter(adjusted_mean_obs, param_number.of.groups)) %>%
    select(-c(qval, sum_mean_obs_counts, num_aggregated_transcripts))
  
  # diagnostic plots
  dispersion <- ggplot(genes.aggregated, aes(x = percent_rank(adjusted_mean_obs), y = -log10(pval))) + 
    geom_point() + 
    ggtitle("dispersion diagram of genes.aggregated with covariate adjusted_mean_obs = sum_mean_obs_counts / num_aggregated_transcripts") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    xlab("percent rank of adjusted_mean_obs") +
    ylab(expression(-log[10](p-value))) 
    
  agg_plot1 <- ggplot(genes.aggregated, aes(x = pval)) +
    geom_histogram(binwidth = 0.025, boundary = 0) +
    xlab("p-values without grouping") +     
    ylab("density")
  
  levels(genes.aggregated$grouping) <- paste0(rep("group ", param_number.of.groups), 1:param_number.of.groups)
  agg_plot2 <- ggplot(genes.aggregated, aes(x = pval)) +
    geom_histogram(binwidth = 0.025, boundary = 0) +
    xlab("p-values") +     
    ylab("density") +
    facet_wrap(~ grouping, nrow = 1) 
  
  plots_agg_mean_obs <- ggarrange(agg_plot1, agg_plot2, nrow = 2, heights = c(3, 1.5))
  
  histograms <- annotate_figure(plots_agg_mean_obs, 
                                 top = text_grob("histograms of genes.aggregated: covariate  = sum_mean_obs_counts / num_aggregated_transcripts",
                                                 color = "black", 
                                                 face = "bold", 
                                                 size = 14))
  
  ggexport(list(dispersion, histograms),
           filename = paste0(output_dir_diagn_plots, "/ihw_genes.aggregated.cov_mean_obs.diagnostic.plots.pdf"), 
           width=14)
  
#######################################################################################################################################
  
# covariates: mean_obs and var_obs in transcripts
transcripts <- read_tsv(".test/results/tables/diffexp/model_X.transcripts.diffexp.tsv")
transcripts_selected <- na.omit(transcripts) %>%
  select(target_id, ens_gene, ext_gene, pval, mean_obs, var_obs) 

transcripts_group_mean <- transcripts_selected %>%
  mutate(grouping = groups_by_filter(mean_obs, param_number.of.groups))

transcripts_group_var <- transcripts_selected %>%
  mutate(grouping = groups_by_filter(var_obs, param_number.of.groups))

# diagnostic plots           
pvalue <- transcripts$pval
cov_mean <- transcripts$mean_obs
cov_var <- transcripts$var_obs

### dispersion plots for transcripts 
dispersion <- rbind(data.frame(pval = pvalue, covariate = percent_rank(cov_mean), title = "mean_obs"),
      data.frame(pval = pvalue, covariate = percent_rank(cov_var), title = "var_obs")) %>%  #data.frame(pval = pvalue, covariate = percent_rank(1/cov_var), title = "var_obs")) %>%
  ggplot(aes(x = covariate, y = -log10(pval))) + 
  geom_point() + 
  facet_grid( . ~ title) +
  ylab(expression(-log[10](p-value)))
  
### histograms of transcripte for mean_obs
transcr_mean_plot1 <- ggplot(transcripts_group_mean, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values without grouping") +     
  ylab("density")

levels(transcripts_group_mean$grouping) <- paste0(rep("group ", param_number.of.groups), 1:param_number.of.groups)
transcr_mean_plot2 <- ggplot(transcripts_group_mean, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values grouped by mean_obs") +     
  ylab("density") +
  facet_wrap(~ grouping, nrow = 1) 

plots_transcr_mean_obs <- ggarrange(transcr_mean_plot1, transcr_mean_plot2, nrow = 2, heights = c(3, 1.5))

histograms_mean_obs <- annotate_figure(plots_transcr_mean_obs, 
                              top = text_grob("histograms of transcripts: covariate  = mean_obs",
                                              color = "black", 
                                              face = "bold", 
                                              size = 14))

### histograms of transcripte for var_obs
transcr_var_plot1 <- ggplot(transcripts_group_var, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values without grouping") +     
  ylab("density")

levels(transcripts_group_var$grouping) <- paste0(rep("group ", param_number.of.groups), 1:param_number.of.groups)
transcr_var_plot2 <- ggplot(transcripts_group_var, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values grouped by var_obs") +     
  ylab("density") +
  facet_wrap(~ grouping, nrow = 1) 

plots_transcr_var_obs <- ggarrange(transcr_var_plot1, transcr_var_plot2, nrow = 2, heights = c(3, 1.5))

histograms_var_obs <- annotate_figure(plots_transcr_var_obs, 
                                       top = text_grob("histograms of transcripts: covariate  = var_obs",
                                                       color = "black", 
                                                       face = "bold", 
                                                       size = 14))

ggexport(list(dispersion, histograms_mean_obs, histograms_var_obs),
         filename = paste0(output_dir_diagn_plots, "/ihw_transcripts.cov_mean_obs_cov_var_obs.diagnostic.plots.pdf"), 
         width=14)

# removing plots and plot specific variables from workspace
rm(transcr_mean_plot1, transcr_mean_plot2, transcr_var_plot1, transcr_var_plot2, plots_agg_mean_obs, plots_transcr_mean_obs,
   plots_transcr_var_obs, agg_plot1, agg_plot2, histograms, histograms_mean_obs, histograms_var_obs, dispersion, transcripts_selected,
   transcripts_group_mean, transcripts_group_var)

#######################################################################################################################################

# ihw implementation

ihw_results_mean <- ihw(pval ~ mean_obs, data = transcripts, alpha = 0.1, nbins = param_number.of.groups)
ihw_mean <- as.data.frame(ihw_results_mean)
ihw_results_var <- ihw(pval ~ var_obs, data = transcripts, alpha = 0.1, nbins = param_number.of.groups)
ihw_var <- as.data.frame(ihw_results_var)

# plot of trends of the covariate 
plot_trend_mean <- plot(ihw_results_mean) + 
  ggtitle("plot of trends of mean_obs as covariate")
plot_trend_var <- plot(ihw_results_var)+ 
  ggtitle("plot of trends of var_obs as covariate")

# plots of decision boundaries for unweighted p-values as a function of the covariate
plot_decision_mean <- plot(ihw_results_mean, what = "decisionboundary") + 
  ggtitle("plot of decision boundaries for unweighted p-values as a function of the covariate mean_obs")
plot_decision_var <- plot(ihw_results_var, what = "decisionboundary") + 
  ggtitle("plot of decision boundaries for unweighted p-values as a function of the covariate var_obs")

# p-values vs. adusted p-values
plot_ihw_pval_mean <- ggplot(ihw_mean, aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + 
  ggtitle("p-values vs. adusted p-values in ihw-analysis of the covariate mean_obs") +
  scale_colour_hue(l = 70, c = 150, drop = FALSE)
  
plot_ihw_pval_var <- ggplot(ihw_var, aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25)  + 
  ggtitle("p-values vs. adusted p-values in ihw-analysis of the covariate mean_obs") +
  scale_colour_hue(l = 70, c = 150, drop = FALSE)

plots_results_mean <- ggarrange(plot_trend_mean, plot_decision_mean, plot_ihw_pval_mean, ncol = 1, nrow = 1, heights = c(2, 1.5))
ggexport(plots_results_mean, filename = str_c(output_dir_result_plots, "/ihw_transcripts.cov_mean_obs.result.plots.pdf"))

plots_results_var <- ggarrange(plot_trend_var, plot_decision_var, plot_ihw_pval_var, ncol = 1, nrow = 1, heights = c(2, 1.5))
ggexport(plots_results_var, filename = str_c(output_dir_result_plots, "/ihw_transcripts.cov_var_obs.result.plots.pdf"))
