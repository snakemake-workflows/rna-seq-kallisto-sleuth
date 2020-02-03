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

so <- sleuth_load(snakemake@input[[1]])

# so <- sleuth_load("results/sleuth/diffexp/model_X.transcripts.diffexp.rds")

#---plot example for testing snakemakerule---
#ggplot(data=iris, aes(x = Sepal.Length, y = Sepal.Width))
#ggsave(snakemake@output[[1]], width = 14)
#--------------------------------------------

pdf(file = snakemake@output[[1]], width = 14)
sleuth::plot_volcano(so, test_type = "wt", which_model = 'full',  
             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
             highlight = NULL)

dev.off()

