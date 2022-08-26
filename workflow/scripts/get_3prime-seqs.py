
import Bio
import json
import re

f = open(snakemake.input["read_length"])
read_length = json.load(f)
f.close()
three_prime_output_file = open(snakemake.output[0], "w")

from Bio import SeqIO
for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
    transcript_location = seq_record.description.split(" ")[2]
    strand = transcript_location.split(":")[5]
    if strand == "1":
        trimmed_seq_postive = seq_record.seq[-read_length:]
        print(">",seq_record.id, sep = "", file = three_prime_output_file)
        print(trimmed_seq_postive, file = three_prime_output_file)
    elif strand == "-1":
        trimmed_seq_negative = seq_record.seq[0:read_length]
        print(">",seq_record.id, sep = "", file = three_prime_output_file)
        print(trimmed_seq_negative, file = three_prime_output_file)

three_prime_output_file.close()