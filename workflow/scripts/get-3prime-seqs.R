# TODO
# 1. filter input fasta (snakemake@input[[1]])
#to canonical transcripts (obtain info from biomart)
# 2. out of the remaining transcripts,
#extract the last n (n <- snakemake@params[["read_length"]])
#bases from the sequence and write into new fasta file (snakemake@output[[1]])
library(biomaRt)
library(dplyr)
library(stringr)
library(rjson)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
    version = snakemake@params[["release"]])
get_canonical_transcripts <-
    getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version",
    "ensembl_transcript_id", "ensembl_transcript_id_version",
    "transcript_is_canonical"), mart = ensembl)
read_length <- fromJSON(file = snakemake@input[["read_length"]])
canonical_ids <- subset(get_canonical_transcripts, transcript_is_canonical == 1)
seq <- getSequence(id = canonical_ids$ensembl_transcript_id_version,
    type = "ensembl_transcript_id_version", seqType = "3utr",
    mart = ensembl)
seq_fil <-
    seq %>% filter(`3utr` != "Sequence unavailable")
seq_fil <- seq_fil %>% mutate(`3utr` = str_remove_all(`3utr`, "AAA"))
seq_fil <- seq_fil %>% mutate(`3utr` = str_remove_all(`3utr`, "TTT"))
#seq_fil <- seq_fil %>% mutate(cdna = str_remove_all(cdna, "AAA"))
seq_fil$`3utr` <- str_sub(seq_fil$`3utr`,
    start = -read_length)
#unava_seq_fil <-
    #seq %>% filter(`3utr` != "Sequence unavailable")
    # if unavailable sequences are to be removed. 
exportFASTA(seq_fil, file = snakemake@output[[1]])