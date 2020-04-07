log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("ggpubr")
library("IHW")

number.of.groups <- 7

gene.data <- na.omit(read_tsv((snakemake@input[[1]])))
level.name <- snakemake@wildcards[["level"]]

# calculate covariate mean_obs for gene.aggregated data
if(level.name == "genes-aggregated") {
  gene.data <- gene.data %>%
    mutate(mean_obs = if_else(num_aggregated_transcripts > 0, (sum_mean_obs_counts / num_aggregated_transcripts), 0))
}

# determine the appropriate number of groups for grouping in ihw calculation for genes_aggregated
tested.number.of.groups <- number.of.groups
ihw_test_grouping <- function(x){
  tryCatch(
    expr = {
      ihw(pval ~ mean_obs, data = gene.data, alpha = 0.1, nbins = x)
      return(TRUE)
    },
    error = function(e) {
      print(str_c("Number of groups was set too hight, trying grouping by ", tested.number.of.groups - 1, " groups."))
      return(FALSE)
    })
}

while(!ihw_test_grouping(tested.number.of.groups) && tested.number.of.groups > 0){
  tested.number.of.groups <- tested.number.of.groups -1
  ihw_test_grouping(tested.number.of.groups)
}

if (tested.number.of.groups < number.of.groups) {
  print(str_c("The chosen number of groups for ", level.name, " dataset was too large, instead IHW was calculated on ", tested.number.of.groups, " groups."))
  number.of.groups <- tested.number.of.groups
}

#select the necessary data
if(level.name == "genes-aggregated") {
  gene.data <- gene.data %>%
    select(-c(qval, sum_mean_obs_counts, num_aggregated_transcripts)) %>%
    mutate(grouping = groups_by_filter(mean_obs, number.of.groups))
} else {
  gene.data <- gene.data %>%
    select(ens_gene, ext_gene, pval, mean_obs) %>%
    mutate(grouping = groups_by_filter(mean_obs, number.of.groups))
}

### diagnostic plots for covariate and grouping
#dispersion plot
dispersion <- ggplot(gene.data, aes(x = percent_rank(mean_obs), y = -log10(pval))) +
  geom_point() +
  ggtitle("IHW: scatter plot of p-values vs. mean of counts as covariate") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  xlab("percent rank of mean_obs") +
  ylab(expression(-log[10](p-value)))

ggsave(snakemake@output[["dispersion"]], dispersion)

#histograms
hist_overall <- ggplot(gene.data, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values without grouping") +
  ylab("density")

levels(gene.data$grouping) <- paste0(rep("group ", number.of.groups), 1:number.of.groups)
hist_groups <- ggplot(gene.data, aes(x = pval)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  xlab("p-values of the individual groups") +
  ylab("density") +
  facet_wrap(~ grouping, nrow = number.of.groups %/% 4 + ifelse((number.of.groups %% 4 != 0), 1, 0))

plots_agg_mean_obs <- ggarrange(hist_overall, hist_groups, nrow = 1)

histograms <- annotate_figure(plots_agg_mean_obs,
                              top = text_grob("IHW: histograms for p-values of mean of counts as covariate",
                                              color = "black",
                                              face = "bold",
                                              size = 12))

ggexport(histograms,
         filename = snakemake@output[["histograms"]],
         width=14)

# ihw calculation
write_tsv(gene.data, snakemake@output[["results"]])
ihw_results_mean <- ihw(pval ~ mean_obs, data = gene.data, alpha = 0.1, nbins = tested.number.of.groups)
ihw_mean <- as.data.frame(ihw_results_mean)
write_tsv(ihw_mean, snakemake@output[["results"]])

### diagnostic plots for ihw calculation
# plot of trends of the covariate
plot_trend_mean <- plot(ihw_results_mean) +
  ggtitle("IHW: Plot of trends of covariate (mean of counts)") +
  theme(plot.title = element_text(size=12))
ggsave(snakemake@output[["trends"]], plot_trend_mean)

# plots of decision boundaries for unweighted p-values as a function of the covariate
plot_decision_mean <- plot(ihw_results_mean, what = "decisionboundary") +
  ggtitle("IHW: Decision boundaries for unweighted p-values vs. mean of counts as covariate") +
  theme(plot.title = element_text(size=12))
ggsave(snakemake@output[["decision"]], plot_decision_mean)

# p-values vs. adusted p-values
plot_ihw_pval_mean <- ggplot(ihw_mean, aes(x = pvalue, y = adj_pvalue, col = group)) +
  geom_point(size = 0.25) +
  ggtitle("IHW: p-values vs. adusted p-values in ihw-analysis of mean of counts as covariate") +
  theme(plot.title = element_text(size=12)) +
  scale_colour_hue(l = 70, c = 150, drop = FALSE)
ggsave(snakemake@output[["adj_pvals"]], plot_ihw_pval_mean)