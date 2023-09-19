log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(biomaRt)
library(tidyverse)

ensembl <- useEnsembl(
    biomart = "genes", dataset = "hsapiens_gene_ensembl",
    version = snakemake@params["release"]
)
transcripts <-
    getBM(
        attributes = c(
            "ensembl_transcript_id_version",
            "transcript_length",
            "strand",
            "transcript_mane_select",
            "chromosome_name"
        ),
        mart = ensembl
    )
cleaned_transcripts <- transcripts %>%
    filter(!str_detect(chromosome_name, "patch|PATCH")) %>%
    select(!chromosome_name) %>%
    rename(transcript = ensembl_transcript_id_version) %>%
    mutate(
        strand = case_match(
            strand,
            1 ~ "+",
            -1 ~ "-"
        ),
        transcript_mane_select = if_else(transcript_mane_select != "", 1, 0)
    )

write_tsv(
    cleaned_transcripts,
    snakemake@output[["all"]]
)

# eventually, use only reads mapping to transcripts from the
# transcript_mane_select set, as these:
# * include one well-supported transcript per protein-coding locus (and only
#   protein-coding transcripts should be targeted by QuantSeq via poly A tail
#   priming)
# * are versioned and largely stable, only allowing updates if absolutely
#   required
# see: https://mart.ensembl.org/info/genome/genebuild/mane.html
write_tsv(
    cleaned_transcripts %>% filter(transcript_mane_select == 1),
    snakemake@output[["mane_select"]]
)
