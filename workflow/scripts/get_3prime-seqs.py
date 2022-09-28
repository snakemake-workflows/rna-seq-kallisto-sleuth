
import Bio
import json
import re
import sys
sys.stderr = open(snakemake.log[0], "w")


three_prime_output_file = open(snakemake.output[0], "w")

from Bio import SeqIO
for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
    transcript_location = seq_record.description.split(" ")[0]
    transcript_id = transcript_location.split("_")[0]
    print(">",transcript_id, sep = "", file = three_prime_output_file)
    print(seq_record.seq, file = three_prime_output_file)
three_prime_output_file.close()