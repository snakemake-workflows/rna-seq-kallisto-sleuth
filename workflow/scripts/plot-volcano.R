log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#---install required packages for local testing without surrounding snakemake-workflow---
#install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("devtools")
#BiocManager::install("pachterlab/sleuth")
#----------------------------------------------------------------------------------------

library("sleuth")
library("tidyverse")

print("Building volcano plot...")



#setwd("/home/tharja/Schreibtisch/RNA-Seq-Project/tmp_rna-seq-kallisto-sleuth/")
 
pval_transcripts <- sleuth_load(snakemake@input[["pvals"]])
conditions <- read.table(snakemake@input[["samples"]])
matrix <- read.table(snakemake@input[["matrix"]], quote="\"", check.names = F, stringsAsFactors=F, header = TRUE) 

# #pval_transcripts <- sleuth_load("results/sleuth/diffexp/model_X.transcripts.diffexp.rds")
# #conditions <- read.table("results/sleuth/samples.tsv")
# #matrix <- read.table("results/tables/tpm-matrix/model_X.tpm-matrix.tsv", quote="\"", check.names = F, stringsAsFactors=F, header = TRUE)

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

  samples <- as_tibble(matrix, header = TRUE)

  samples_condition_a <- samples %>%
    select(condition_a) 

  samples_condition_b <- samples %>%
    select(condition_b)

#}
 

 
  fold_change <- tibble(target_id=samples$target_id, ext_gene=samples$ext_gene, 
                        FC=rowSums(samples_condition_a)/rowSums(samples_condition_b),
                        FC_reverse=rowSums(samples_condition_b)/rowSums(samples_condition_a))
  p_vals <- pval_transcripts %>%
    select(target_id, ext_gene, pval, qval)
    
  log2FC <- fold_change %>%
    select(target_id, ext_gene, log2FC=FC) %>%
    mutate(log2FC=log2(log2FC)) 
  
  plot_data <- full_join(log2FC, p_vals, by=c("target_id", "ext_gene")) %>%
    mutate(pval=-log10(pval))
  
  # plot_data <- plot_data[is.finite(plot_data$log2FC) && is.finite(plot_data$pval) && is.finite(plot_data$qval)]

  # print(plot_data$log2FC)
  # print(plot_data$pval)
  #plot_data <- na.omit(plot_data)
  # plot_data <- plot_data%>%mutate(sapply(c(log2FC,pval,qval),rm_nan_inf(x)))
  # rm_nan_inf <- function(x){
  #   plot_data$x[!(is.finite(plot_data$x))] <- 0
  #   plot_data$x[is.na(plot_data$x)] <- 0
  #   plot_data$x[is.nan(plot_data$x)] <- 0
  # }
  
  # plot_data <- plot_data[is.finite(plot_data$log2FC)]
 plot_data$log2FC[!(is.finite(plot_data$log2FC))] <- 0
 plot_data$log2FC[is.na(plot_data$log2FC)] <- 0
 plot_data$log2FC[is.nan(plot_data$log2FC)] <- 0

  # plot_data <- plot_data[is.finite(plot_data$pval)]
 plot_data$pval[!(is.finite(plot_data$pval))] <- 0
 plot_data$pval[is.na(plot_data$pval)] <- 0
 plot_data$pval[is.nan(plot_data$pval)] <- 0

  # plot_data <- plot_data[is.finite(plot_data$qval)]
 plot_data$qval[!(is.finite(plot_data$qval))] <- 0
 plot_data$qval[is.na(plot_data$qval)] <- 0
 plot_data$qval[is.nan(plot_data$qval)] <- 0
  
  # plot_data$log2FC[which(is.nan(plot_data$log2FC))] = Inf
  # plot_data$log2FC[which(is.na(plot_data$log2FC))] = 0
  # 
  # plot_data$pval[which(is.nan(plot_data$pval))] = Inf
  # plot_data$pval[which(is.na(plot_data$pval))] = 0
  # 
  # plot_data$qval[which(is.nan(plot_data$qval))] = Inf
  # plot_data$qval[which(is.na(plot_data$qval))] = 0
  
  # cond1 <- df$sub == 1 & df$day == 2
  # 
  # cond2 <- df$sub == 3 & df$day == 4
  # 
  # df <- df[!(cond1 | cond2),]
  # 
  # is.na(plot_data$log2FC)
  

  #plot_data <- drop_na(plot_data)
  
  plot_data <- plot_data %>%
    mutate(colour=factor(case_when(qval<.05 ~"most significant",
                            abs(log2FC)>1 ~"significant",
                            (qval<.05 && abs(log2FC)>1) ~"low significant",
                            (qval>.05 || abs(log2FC)<1)~"not significant")))
  interc_x <- c(-2,2) 
  interc_y <- c(0.5,1)
  max_x <- max(abs(plot_data$log2FC[is.finite(plot_data$log2FC)]), na.rm=T)
  max_y <- max(abs(plot_data$pval[is.finite(plot_data$pval)]), na.rm = T)
  if(is.numeric(snakemake@params[["range_log2FC"]])){
    max_x <- snakemake@params[["range_log2FC"]]
    interc_x <- c(-max_x/4, max_x/4)
  }
  
  if(is.numeric(snakemake@params[["range_pval"]])){
    max_y <- snakemake@params[["range_pval"]]
    interc_y <- c(0, max_y/4)
  }
   # print(max_x)
   # print(max_y)
   
   # print(plot_data$log2FC)
   # print(plot_data$pval)

  #---plot example for testing snakemakerule---
  #pdf(file = snakemake@output[[1]], width = 14)
  #example <- c(1, 2, 3, 7, 4, 9, 5)
  #plot(example)
  #dev.off()
  #--------------------------------------------
  write_tsv(fold_change, snakemake@output[["foldchange"]])
  pdf(file = snakemake@output[["volcanoplot"]], width = 14)
  #par()

  #pdf(file = "results/plots/volcano/model_X.volcano.png", width = 14)
    volcano <- ggplot(plot_data, aes(x=plot_data$log2FC, y=plot_data$pval))+
      geom_point(aes(color=plot_data$colour)) +
      ggtitle("Volcano-Plot")+
      xlab("log2 fold change")+
      ylab("p-values as -log10 adjusted")+
      xlim(-max_x, max_x)+
      ylim(0, max_y)+
      scale_color_manual(name="significance",
                         values = c("most significant"="red",
                                    "significant"="orange",
                                    "low significant"="green",
                                    "not significant"="black"))+
      geom_vline(xintercept = interc_x, color="blue", alpha=0.5)+
      geom_hline(yintercept = interc_y, color="red", alpha=0.5)+
      geom_vline(xintercept = 0, color="red", alpha=1)
  print(volcano)
  #---plot example for testing snakemakerule---
  #pdf(file = snakemake@output[[1]], width = 14)
  #example <- c(1, 2, 3, 7, 4, 9, 5)
  #plot(example)
  #dev.off()
  #--------------------------------------------
  # ggsave(filename=snakemake@output[[1]], width = 14)
# #  write_tsv("foldchange", snakemake@output[["foldchange"]])
   # print(volcano)
  dev.off()
  
  #ggsave(filename=snakemake@output[[1]], plot = last_plot(), device="pdf", width = 14)
}
  

#########################################################################
#########################################################################

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
