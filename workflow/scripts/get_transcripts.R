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
            "transcript_mane_select",
            "chromosome_name",
            "transcript_length",
            "strand"
        ),
        mart = ensembl
    )
cleaned_transcripts <- transcripts %>%
    filter(!str_detect(chromosome_name, "patch|PATCH")) %>%
    mutate(
        strand = case_match(
            1 ~ "+",
            -1 ~ "-"
        ),
        transcript_mane_select = if_else(transcript_mane_select != "", 1, 0)
    )
    # add empty columns necessary for valid BED file and sort accordingly, see:
    # https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    add_column(
        start = 0,
        name = ""
    ) %>%
    # here, we hijack the `score` field of the BED format to encode the MANE membership
    select(
        ensembl_transcript_id_version,
        start,
        transcript_length,
        name,
        transcript_mane_select,
        strand
    )

write_tsv(
    cleaned_transcripts,
    snakemake@output[["all"]],
    col_names = FALSE
)

# use only those annotated as transcript_mane_select, as these:
# * include one well-supported transcript per protein-coding locus (and only protein-coding transcripts should be targeted by QuantSeq via poly A tail priming)
# * are versioned and largely stable, only allowing updates if absolutely required
# see: https://mart.ensembl.org/info/genome/genebuild/mane.html
write_tsv(
    cleaned_transcripts %>% filter(transcript_mane_select == 1),
    snakemake@output[["mane_select"]],
    col_names = FALSE
)
