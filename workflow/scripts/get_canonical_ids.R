log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(biomaRt)
library(dplyr)
library(tidyverse)

ensembl <- useEnsembl(
    biomart = "genes", dataset = "hsapiens_gene_ensembl",
    version = snakemake@params["release"]
)
get_transcripts_ids <-
    getBM(
        attributes = c(
            "ensembl_transcript_id_version",
            "transcript_is_canonical", "transcript_mane_select", "chromosome_name"
        ),
        mart = ensembl
    )
# transcript_mane_select- manually curated primary transcript as it has been encoded in NCBI
canonical_ids <- get_transcripts_ids %>%
    select(
        ensembl_transcript_id_version, transcript_is_canonical,
        transcript_mane_select, chromosome_name
    ) %>%
    filter(!str_detect(chromosome_name, "patch|PATCH")) %>%
    filter(str_detect(transcript_mane_select, "")) %>%
    subset(transcript_is_canonical == 1)
write.csv(canonical_ids$ensembl_transcript_id_version,
    file = snakemake@output[[1]], quote = FALSE, row.names = FALSE
)
