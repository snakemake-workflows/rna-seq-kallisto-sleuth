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

print("building volcano plot...")

so <- sleuth_load(snakemake@input[["so"]])
pdf(file = snakemake@output[[1]])

example <- c(1, 2, 3, 7, 4, 9, 5)

# Graph the cars vector with all defaults
plot(example)
#plot_volcano(so, test, test_type = "wt", which_model = snakemake@wildcards[["model"]],
#             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
#             highlight = NULL)

dev.off()

#so <- sleuth_load("results/sleuth/all.rds")
#model <- sleuth_load("results/sleuth/model_X.rds")
#model <- sleuth_fit(so, as.formula(model[["full"]]), 'full')
#plot_volcano(so, test, test_type = "wt", which_model = 'full',
#             sig_level = 0.1, point_alpha = 0.2, sig_color = "red",
#             highlight = NULL)