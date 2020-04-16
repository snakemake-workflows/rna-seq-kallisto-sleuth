log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("ggpubr")
library("IHW")

number_of_groups <- 7

gene_data <- read_tsv(snakemake@input[[1]])
level_name <- snakemake@wildcards[["level"]]

# calculate covariate mean_obs for gene.aggregated data
if(level_name == "genes-aggregated") {
  gene_data <- gene_data %>%
    mutate(mean_obs = if_else(num_aggregated_transcripts > 0, sum_mean_obs_counts / num_aggregated_transcripts, 0), ens_gene = target_id)
}

# determine the appropriate number of groups for grouping in ihw calculation
tested_number_of_groups <- number_of_groups
ihw_test_grouping <- function(x){
  tryCatch(
    expr = {
      ihw(pval ~ mean_obs, data = gene_data, alpha = 0.1, nbins = x)
      return(TRUE)
    },
    error = function(e) {
      print(str_c("Number of groups was set too hight, trying grouping by ", tested_number_of_groups - 1, " groups."))
      return(FALSE)
    })
}

while(!ihw_test_grouping(tested_number_of_groups) && tested_number_of_groups > 0){
  tested_number_of_groups <- tested_number_of_groups -1
  ihw_test_grouping(tested_number_of_groups)
}

if (tested_number_of_groups < number_of_groups) {
  print(str_c("The chosen number of groups for ", level_name, " dataset was too large, instead IHW was calculated on ", tested_number_of_groups, " groups."))
  number_of_groups <- tested_number_of_groups
}

gene_data <- gene_data %>%
  drop_na(pval, mean_obs) %>%
  select(ens_gene, ext_gene, pval, mean_obs) %>%
  mutate(grouping = groups_by_filter(mean_obs, number_of_groups))

### diagnostic plots for covariate and grouping
#dispersion plot
dispersion <- ggplot(gene_data, aes(x = percent_rank(mean_obs), y = -log10(pval))) +
  geom_point() +
  ggtitle("IHW: scatter plot of p-values vs. mean of counts") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  xlab("percent rank of mean_obs") +
  ylab(expression(-log[10](p-value)))

ggsave(snakemake@output[["dispersion"]], dispersion)

#histograms
hist_overall <- ggplot(gene_data, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values without grouping") +
  ylab("density")

hist_groups <- ggplot(gene_data, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values of the individual groups") +
  ylab("density") +
  facet_wrap(~ grouping, nrow = ceiling(number_of_groups / 4)) 

plots_agg_mean_obs <- ggarrange(hist_overall, hist_groups, nrow = 1)

histograms <- annotate_figure(plots_agg_mean_obs,
                              top = text_grob("IHW: histograms for p-values of mean of counts",
                                              color = "black",
                                              face = "bold",
                                              size = 12))

ggexport(histograms,
         filename = snakemake@output[["histograms"]],
         width=14)

# ihw calculation
ihw_results_mean <- ihw(pval ~ mean_obs, data = gene_data, alpha = 0.1, nbins = tested_number_of_groups)

# merging ens_gene-IDs and ext_gene-names
ihw_mean <- as.data.frame(ihw_results_mean) %>% 
  # TODO remove ugly hack if ihw in future allows annotation columns (feature requested here: https://support.bioconductor.org/p/129972/)
  right_join(gene_data, by = c(pvalue = "pval", covariate = "mean_obs", group = "grouping")) %>%
  unique() %>%
  select(ens_gene, ext_gene, everything())

write_tsv(ihw_mean, snakemake@output[["results"]])

### diagnostic plots for ihw calculation
# plot of trends of the covariate
plot_trend_mean <- plot(ihw_results_mean) +
  ggtitle("IHW: Plot of trends of mean of counts") +
  theme(plot.title = element_text(size=12))
ggsave(snakemake@output[["trends"]], plot_trend_mean)

# plots of decision boundaries for unweighted p-values as a function of the covariate
plot_decision_mean <- plot(ihw_results_mean, what = "decisionboundary") +
  ggtitle("IHW: Decision boundaries for unweighted p-values vs. mean of counts") +
  theme(plot.title = element_text(size=12))
ggsave(snakemake@output[["decision"]], plot_decision_mean)

# p-values vs. adusted p-values
plot_ihw_pval_mean <- ggplot(ihw_mean, aes(x = pvalue, y = adj_pvalue, col = group)) +
  geom_point(size = 0.25) +
  ggtitle("IHW: raw p-values vs. adusted p-values from ihw-analysis") +
  theme(plot.title = element_text(size=12)) +
  scale_colour_hue(l = 70, c = 150, drop = FALSE)
ggsave(snakemake@output[["adj_pvals"]], plot_ihw_pval_mean)
