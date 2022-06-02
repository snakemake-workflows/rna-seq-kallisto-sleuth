# TODO
# 1. filter input fasta (snakemake@input[[1]]) to canonical transcripts (obtain info from biomart)
# 2. out of the remaining transcripts, extract the last n (n <- snakemake@params[["read_length"]]) bases from the sequence and write into new fasta file (snakemake@output[[1]])