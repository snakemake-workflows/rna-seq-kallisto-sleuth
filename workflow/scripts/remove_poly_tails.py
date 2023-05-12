import Bio
import sys
import re
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as transcript_clean_cdna_fasta:
    for seq_record in SeqIO.parse(snakemake.input["ref_fasta"], "fasta"):
        transcript_location = seq_record.description.split(" ")[2]
        strand = transcript_location.split(":")[5]
        if strand == "1":
            polyrem_seq = re.sub("TTTT+$|AAAA+$", "", str(seq_record.seq))
            print(
                ">",
                seq_record.id,
                "_1 ",
                seq_record.description,
                sep="",
                file=transcript_clean_cdna_fasta,
            )
            print(polyrem_seq, file=transcript_clean_cdna_fasta)
        elif strand == "-1":
            polyrem_seq = re.sub("^TTTT+|^AAAA+", "", str(seq_record.seq))
            print(
                ">",
                seq_record.id,
                "_-1 ",
                seq_record.description,
                sep="",
                file=transcript_clean_cdna_fasta,
            )
            print(polyrem_seq, file=transcript_clean_cdna_fasta)
