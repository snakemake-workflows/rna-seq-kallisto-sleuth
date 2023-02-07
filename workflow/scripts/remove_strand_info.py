import Bio
import sys
import re
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as three_prime_output_file:
    for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
        transcript_location = seq_record.description.split(" ")[0]
        # Split the transcript id by `_` to filter the strand info from transcript id
        transcript_id = transcript_location.split("_")[0]
        print(">", transcript_id, sep="", file=three_prime_output_file)
        print(seq_record.seq, file=three_prime_output_file)
