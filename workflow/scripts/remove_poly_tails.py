import Bio
import json
import re

transcript_clean_cdna_fasta = open(snakemake.output[0], "w")

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
        #trimmed_seq_postive = polyrem_seq[-5:]
        #print(">",seq_record.id, sep = "")
        #print(trimmed_seq_postive)
        print(">",seq_record.description, sep = "", file = transcript_clean_cdna_fasta)
        print(polyrem_seq, file = transcript_clean_cdna_fasta)
    elif strand == "-1":
        polyrem_seq = re.sub('^T+|^A+', '', str(seq_record.seq))
        #trimmed_seq_negative = polyrem_seq[0:5]
        #print(">",seq_record.id, sep = "")
        #print(trimmed_seq_negative)
        print(">",seq_record.description, sep = "", file = transcript_clean_cdna_fasta)
        print(polyrem_seq, file = transcript_clean_cdna_fasta)
transcript_clean_cdna_fasta.close()