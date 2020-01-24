#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

#---install required packages for local testing without surrounding snakemake-workflow---
#install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("devtools")
#BiocManager::install("pachterlab/sleuth")
#----------------------------------------------------------------------------------------

library("sleuth")
library("tidyverse")

print("Building volcano plot...")

setwd("/home/tharja/Schreibtisch/RNA-Seq-Project/tmp_rna-seq-kallisto-sleuth/")

#pval_transcripts <- sleuth_load(snakemake@input[["pvals"]])
#conditions <- read.table(snakemake@input[["samples"]])
#matrix <- sleuth_load(snakemake@input[["matrix"]]) 

pval_transcripts <- sleuth_load("results/sleuth/diffexp/model_X.transcripts.diffexp.rds")
conditions <- read.table("results/sleuth/samples.tsv")
matrix <- read.table("results/tables/tpm-matrix/model_X.tpm-matrix.tsv", quote="\"", check.names = F, stringsAsFactors=F, header = TRUE)

colnames(matrix)[1] = "target_id"
matrix <- as_tibble(matrix)

levels_condition <- levels(conditions$condition)

if(length(levels_condition)!=2) {
  print("Building volcano plot failed! There are too much different conditions!")
} else {
  condition_a_name <- levels_condition[1]
  condition_b_name <- levels_condition[2]

  condition_a <- as.character(factor(conditions$sample[conditions$condition == condition_a_name]))
  condition_b <- as.character(factor(conditions$sample[conditions$condition == condition_b_name]))
}
  samples <- as_tibble(matrix, header = TRUE)

  samples_condition_a <- samples %>%
    select(condition_a) 

  samples_condition_b <- samples %>%
    select(condition_b)

  
  fold_change <- tibble(target_id=samples$target_id, ext_gene=samples$ext_gene, 
                        FC=rowSums(samples_condition_a)/rowSums(samples_condition_b),
                        FC_reverse=rowSums(samples_condition_b)/rowSums(samples_condition_a))
  p_vals <- pval_transcripts %>%
    select(target_id, ext_gene, pval, qval)
    
  log2FC <- fold_change %>%
    select(target_id, ext_gene, log2FC=FC) %>%
    mutate(log2FC=log2(log2FC)) 
  
  plot_data <- full_join(log2FC, p_vals, by=c("target_id", "ext_gene")) %>%
    mutate(pval=-log10(pval)) %>%
    filter(is.finite(log2FC)) %>%
    filter(is.finite(pval)) %>%
    filter(is.finite(qval))

  plot_data <- na.omit(plot_data)
  
  plot_data <- plot_data %>%
    mutate(colour=factor(case_when(qval<.05 ~"most significant",
                            abs(log2FC)>1 ~"significant",
                            (qval<.05 && abs(log2FC)>1) ~"low significant",
                            (qval>.05 || abs(log2FC)<1)~"not significant")))
  
#  pdf(file = snakemake@output[[1]])
    ggplot(plot_data, aes(x=plot_data$log2FC, y=plot_data$pval))+
      geom_point(aes(color=plot_data$colour))+
      ggtitle("Volcano-Plot")+
      xlab("log2 fold change")+
      ylab("p-values (-log10 adjusted)")+
      xlim(-max(abs(plot_data$log2FC)), max(abs(plot_data$log2FC)))+
      ylim(0, max(abs(plot_data$pval)))+
      scale_color_manual(name="significance",
                         values = c("most significant"="red",
                                    "significant"="orange",
                                    "low significant"="green",
                                    "not significant"="black"))+
      geom_vline(xintercept = c(-2,2), color="blue", alpha=0.5)+
      geom_hline(yintercept = c(0.1, 0.5), color="red", alpha=0.5)
    

                        
#  write_results("foldchange", snakemake@output[["foldchange"]])
  
#  dev.off()
#}

#ltr_data <- sleuth_load(snakemake@input[["ltr"]])
#wt_data <- sleuth_load(snakemake@input[["wt"]])
#smatrix <- sleuth_load(snakemake@input[["matrix"]])

#---for testing---
#so <- sleuth_load("results/sleuth/all.rds")
#model <- sleuth_load("results/sleuth/model_X.rds")
#model <- sleuth_fit(so, as.formula(model[["full"]]), 'full')
#plot_volcano(so, test, test_type = "wt", which_model = 'full',
#             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
#             highlight = NULL)


#so <- sleuth_load("results/sleuth/all.rds")
#transcripts <- sleuth_load("results/sleuth/diffexp/model_X.transcripts.diffexp.rds")
#aggregated <- sleuth_load("results/sleuth/diffexp/model_X.genes-aggregated.diffexp.rds")
#most_sig <- sleuth_load("results/sleuth/diffexp/model_X.genes-mostsigtrans.diffexp.rds")
#model_rds <-sleuth_load("results/sleuth/model_X.rds")

#View(all_rds)
#View(transcripts)
#View(aggregated)
#View(most_sig)
#View(model_rds)

#so <- model_rds

#so <- sleuth_fit(so, as.formula(model_X[["full"]]), 'full')
#so <- sleuth_fit(so, as.formula(model[["reduced"]]), 'reduced')
#sleuth_results(so, covariate, "wt", show_all = TRUE, pval_aggregate = FALSE)
#sleuth_wt(so)
#-----------------

#sleuth::plot_volcano(so, which_model = 'full', test_type = "wt", 
#             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
#             highlight = NULL)

#pdf(file = snakemake@output[[1]])

#---plot example for testing snakemakerule---
#example <- c(1, 2, 3, 7, 4, 9, 5)
#plot(example)
#--------------------------------------------

#plot_volcano(so, test, test_type = "wt", which_model = snakemake@wildcards[["model"]],
#             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
#             highlight = NULL)

#all <- sleuth_load(snakemake@input[[1]])
#so <- sleuth_load("results/sleuth/all.rds")
