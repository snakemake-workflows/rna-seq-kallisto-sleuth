log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("sleuth")

so <- sleuth_load(snakemake@input[["so"]])

diffexp <- read.table(snakemake@input[["diffexp"]], sep = "\t", header = TRUE)

pdf(file = snakemake@output[[1]], width = 14)
plot_transcript_heatmap(so, transcripts = diffexp$target_id[1:20]) 
# TODO, once pachterlab/sleuth#214 is merged, add this to get gene names
# , labels_row = diffexp$ext_gene[1:20])
dev.off()
