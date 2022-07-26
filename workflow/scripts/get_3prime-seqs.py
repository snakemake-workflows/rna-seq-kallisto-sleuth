
import Bio
import json
f = open(snakemake.input["read_length"])
read_length = json.load(f)
f.close()
three_prime_output_file = open(snakemake.output[0], "w")
from Bio import SeqIO
for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
#for seq_record in SeqIO.parse("test1.fa", "fasta"):
    header_info = seq_record.description.split(":")
    strand = header_info[5].replace(" gene","")
    strand.replace(" ","")
    #print("****"+strand+"***")
    if strand == "1":
        seqs = seq_record.seq.rstrip("A|T")
        trimmed_seq_postive = seqs[-read_length:]
        print(">",seq_record.id, sep = "", file = three_prime_output_file)
        print(trimmed_seq_postive, file = three_prime_output_file)
    elif strand == "-1":
        seqs = seq_record.seq.lstrip("A|T")
        trimmed_seq_negative = seqs[0:read_length]
        print(">",seq_record.id, sep = "", file = three_prime_output_file)
        print(trimmed_seq_negative, file = three_prime_output_file)
three_prime_output_file.close()