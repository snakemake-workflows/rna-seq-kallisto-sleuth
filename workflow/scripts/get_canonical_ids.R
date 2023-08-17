log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(biomaRt)
library(tidyverse)

ensembl <- useEnsembl(
    biomart = "genes", dataset = "hsapiens_gene_ensembl",
    version = snakemake@params["release"]
)
get_transcripts_ids <-
    getBM(
        attributes = c(
            "ensembl_transcript_id_version",
            "transcript_is_canonical",
            "transcript_mane_select",
            "chromosome_name",
            "transcript_length",
            "strand"
        ),
        mart = ensembl
    )
# transcript_mane_select- manually curated primary transcript as it has been encoded in NCBI
canonical_ids <- get_transcripts_ids %>%
    filter(!str_detect(chromosome_name, "patch|PATCH")) %>%
    filter(str_detect(transcript_mane_select, "")) %>%
    filter(transcript_is_canonical == 1) %>%
    # add columns necessary for valid BED file and sort accordingly, see:
    # https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    add_column(
        start = 0
    ) %>%
    # Here, we hijack the `name` field of the BED format to encode the strand,
    # as biomart gives the strand as `1` vs `-1`, as opposed to BED's `+` vs. `-`.
    select(
        ensembl_transcript_id_version,
        start,
        transcript_length,
        strand
    )

write_tsv(
    canonical_ids,
    snakemake@output[[1]],
    col_names = FALSE
)
