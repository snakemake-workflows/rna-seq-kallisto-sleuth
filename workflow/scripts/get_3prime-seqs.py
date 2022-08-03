
import Bio
import json
import re

f = open(snakemake.input["read_length"])
read_length = json.load(f)
f.close()
three_prime_output_file = open(snakemake.output[0], "w")

from Bio import SeqIO
for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
#for seq_record in SeqIO.parse("test1.fa", "fasta"):
    transcript_location = seq_record.description.split(" ")[2]
    #print(seq_record.seq)
    #polyrem_seq = re.sub('T+.$|A+.$', '', str(seq_record.seq))
    #print(polyrem_seq)
    strand = transcript_location.split(":")[5]
    if strand == "1":
        polyrem_seq = re.sub('T+$|A+$', '', str(seq_record.seq))
        trimmed_seq_postive = polyrem_seq[-read_length:]
        #trimmed_seq_postive = polyrem_seq[-5:]
        #print(">",seq_record.id, sep = "")
        #print(trimmed_seq_postive)
        print(">",seq_record.id, sep = "", file = three_prime_output_file)
        print(trimmed_seq_postive, file = three_prime_output_file)
    elif strand == "-1":
        polyrem_seq = re.sub('^T+|^A+', '', str(seq_record.seq))
        trimmed_seq_negative = polyrem_seq[0:read_length]
        #trimmed_seq_negative = polyrem_seq[0:5]
        #print(">",seq_record.id, sep = "")
        #print(trimmed_seq_negative)
        print(">",seq_record.id, sep = "", file = three_prime_output_file)
        print(trimmed_seq_negative, file = three_prime_output_file)
three_prime_output_file.close()