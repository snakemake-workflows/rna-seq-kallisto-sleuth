library(biomaRt)
library(tximport)
library(pheatmap)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
    version = snakemake@params[["release"]])
get_transcripts_ids <-
    getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version",
                "ensembl_transcript_id", "ensembl_transcript_id_version"),
                    mart = ensembl)
subset_trans_ids <-
    data.frame(get_transcripts_ids$ensembl_transcript_id_version,
        get_transcripts_ids$ensembl_gene_id)
kallisto_files <- file.path(snakemake@input[["kallisto_path"]], "abundance.h5")
names(kallisto_files) <- list.dirs(snakemake@params[["sample_names"]],
    recursive = FALSE, full.names = FALSE)
txi_kallisto_g <-
    tximport(kallisto_files, type = "kallisto", tx2gene = subset_trans_ids)

print(head(txi_kallisto_g$abundance))
kallsito_file <- data.frame(txi_kallisto_g$abundance)
write.csv(txi_kallisto_g$abundance, snakemake@output[["matrix_file"]],
    quote = FALSE)
vargenes <-
    apply(txi_kallisto_g$abundance, 1, var)
selectedgenes <-
    names(vargenes[order(vargenes, decreasing = TRUE)][1:50])
png(snakemake@output[["png"]], w = 908, h = 953, pointsize = 200)
pheatmap(txi_kallisto_g$abundance[selectedgenes, ], scale = "row")
dev.off()
